import pytest
import os
import pandas as pd
import tarfile
import hashlib
from pathlib import Path
from tempfile import NamedTemporaryFile
from unittest.mock import patch, mock_open
from platformdirs import user_cache_dir
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from plannotate import resources as rsc


class TestResourcePaths:
    """Tests for resource path functions."""

    def test_get_resource(self):
        """Test get_resource function."""
        result = rsc.get_resource("data", "test.yml")
        assert "data/test.yml" in result
        assert result.endswith("data/test.yml")

    def test_get_yaml_path(self):
        """Test get_yaml_path function."""
        result = rsc.get_yaml_path()
        assert "databases.yml" in result

    def test_get_image(self):
        """Test get_image function."""
        result = rsc.get_image("test.png")
        assert "images/test.png" in result

    def test_get_template(self):
        """Test get_template function."""
        result = rsc.get_template("test.html")
        assert "templates/test.html" in result

    def test_get_example_fastas(self):
        """Test get_example_fastas function."""
        result = rsc.get_example_fastas()
        assert "fastas" in result

    def test_get_details(self):
        """Test get_details function."""
        result = rsc.get_details("test.csv")
        assert "test.csv" in result


class TestFileValidation:
    """Tests for file validation functions."""

    def test_get_name_ext(self):
        """Test get_name_ext function."""
        name, ext = rsc.get_name_ext("/path/to/test.fasta")
        assert name == "test"
        assert ext == ".fasta"

        name, ext = rsc.get_name_ext("simple.gbk")
        assert name == "simple"
        assert ext == ".gbk"

    def test_validate_sequence_valid(self):
        """Test validate_sequence with valid DNA."""
        valid_seq = "ATCGATCGATCG"
        rsc.validate_sequence(valid_seq)  # Should not raise

    def test_validate_sequence_with_iupac(self):
        """Test validate_sequence with IUPAC codes."""
        iupac_seq = "ATCGATCGRYWSMKHBVDN"
        rsc.validate_sequence(iupac_seq)  # Should not raise

    def test_validate_sequence_invalid_chars(self):
        """Test validate_sequence with invalid characters."""
        invalid_seq = "ATCGATCGZXQ"
        with pytest.raises(ValueError, match="invalid characters"):
            rsc.validate_sequence(invalid_seq)

    def test_validate_sequence_too_long(self):
        """Test validate_sequence with sequence too long."""
        long_seq = "A" * (rsc.MAX_PLAS_SIZE + 1)
        with pytest.raises(ValueError, match="too large"):
            rsc.validate_sequence(long_seq, max_length=rsc.MAX_PLAS_SIZE)

    def test_validate_file_fasta_valid(self):
        """Test validate_file with valid FASTA."""
        fasta_content = ">test\nATCGATCGATCG\n"
        
        with NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
            tmp.write(fasta_content)
            tmp.flush()
            
            try:
                result = rsc.validate_file(tmp.name, '.fasta')
                assert result == "ATCGATCGATCG"
            finally:
                os.unlink(tmp.name)

    def test_validate_file_fasta_multiple_entries(self):
        """Test validate_file with FASTA containing multiple entries."""
        fasta_content = ">test1\nATCG\n>test2\nGCTA\n"
        
        with NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
            tmp.write(fasta_content)
            tmp.flush()
            
            try:
                with pytest.raises(ValueError, match="many entries"):
                    rsc.validate_file(tmp.name, '.fasta')
            finally:
                os.unlink(tmp.name)

    def test_validate_file_fasta_malformed(self):
        """Test validate_file with malformed FASTA."""
        malformed_content = "not a fasta file"
        
        with NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
            tmp.write(malformed_content)
            tmp.flush()
            
            try:
                with pytest.raises(ValueError, match="Malformed fasta"):
                    rsc.validate_file(tmp.name, '.fasta')
            finally:
                os.unlink(tmp.name)

    def test_validate_file_genbank_valid(self):
        """Test validate_file with valid GenBank."""
        # Create a minimal GenBank record
        record = SeqRecord(Seq("ATCGATCGATCG"), id="test", description="test sequence")
        
        with NamedTemporaryFile(mode='w', suffix='.gbk', delete=False) as tmp:
            SeqIO.write(record, tmp.name, 'genbank')
            
            try:
                result = rsc.validate_file(tmp.name, '.gbk')
                assert result == "ATCGATCGATCG"
            finally:
                os.unlink(tmp.name)

    def test_validate_file_invalid_extension(self):
        """Test validate_file with invalid file extension."""
        with NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
            tmp.write("test")
            tmp.flush()
            
            try:
                with pytest.raises(ValueError, match="must be a FASTA or GenBank"):
                    rsc.validate_file(tmp.name, '.txt')
            finally:
                os.unlink(tmp.name)


class TestYamlParsing:
    """Tests for YAML configuration parsing."""

    def test_get_yaml(self):
        """Test get_yaml function."""
        yaml_content = """
database1:
  method: diamond
  location: /path/to/db
  priority: 1
  parameters:
    - "--sensitive"
    - "--max-target-seqs"
    - "500"
database2:
  method: infernal
  location: /path/to/rfam
  priority: 2
"""
        
        with NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as tmp:
            tmp.write(yaml_content)
            tmp.flush()
            
            try:
                result = rsc.get_yaml(tmp.name)
                
                assert 'database1' in result
                assert 'database2' in result
                assert result['database1']['method'] == 'diamond'
                assert result['database1']['db_loc'] == '/path/to/db'
                assert result['database1']['parameters'] == '--sensitive --max-target-seqs 500'
                assert result['database2']['parameters'] == ''  # No parameters
            finally:
                os.unlink(tmp.name)


class TestGenBankGeneration:
    """Tests for GenBank file generation."""

    def test_get_seq_record_empty_df(self):
        """Test get_seq_record with empty DataFrame."""
        empty_df = pd.DataFrame()
        sequence = "ATCGATCGATCG"
        
        record = rsc.get_seq_record(empty_df, sequence, is_linear=True)
        
        assert str(record.seq) == sequence
        assert record.annotations['topology'] == 'linear'
        assert 'pLannotate' in record.annotations['comment']

    def test_get_seq_record_with_features(self):
        """Test get_seq_record with features."""
        df = pd.DataFrame({
            'qstart': [0, 5],
            'qend': [4, 10],
            'qlen': [15, 15],
            'sframe': [1, -1],
            'Feature': ['gene1', 'gene2'],
            'Type': ['CDS', 'misc_feature'],
            'db': ['test_db', 'test_db'],
            'pident': [95.0, 88.0],
            'percmatch': [80.0, 75.0],
            'fragment': [False, True]
        })
        sequence = "ATCGATCGATCGATC"
        
        record = rsc.get_seq_record(df, sequence, is_linear=False)
        
        assert str(record.seq) == sequence
        assert record.annotations['topology'] == 'circular'
        assert len(record.features) == 2
        
        # Check feature properties
        feature1 = record.features[0]
        assert feature1.qualifiers['label'][0] == 'gene1'
        assert feature1.type == 'CDS'
        
        feature2 = record.features[1]
        assert feature2.qualifiers['label'][0] == 'gene2 (fragment)'
        assert feature2.type == 'misc_feature'

    def test_get_seq_record_origin_crossing(self):
        """Test get_seq_record with origin-crossing feature."""
        df = pd.DataFrame({
            'qstart': [10],
            'qend': [5],  # End before start indicates origin crossing
            'qlen': [15],
            'sframe': [1],
            'Feature': ['origin_gene'],
            'Type': ['CDS'],
            'db': ['test_db'],
            'pident': [95.0],
            'percmatch': [80.0],
            'fragment': [False]
        })
        sequence = "ATCGATCGATCGATC"
        
        record = rsc.get_seq_record(df, sequence)
        
        assert len(record.features) == 1
        feature = record.features[0]
        # Should create a compound location for origin crossing
        assert hasattr(feature.location, 'parts') or feature.location.start > feature.location.end

    def test_get_gbk(self):
        """Test get_gbk function."""
        df = pd.DataFrame({
            'qstart': [0],
            'qend': [10],
            'qlen': [15],
            'sframe': [1],
            'Feature': ['test_gene'],
            'Type': ['CDS'],
            'db': ['test_db'],
            'pident': [95.0],
            'percmatch': [80.0],
            'fragment': [False]
        })
        sequence = "ATCGATCGATCGATC"
        
        gbk_text = rsc.get_gbk(df, sequence, is_linear=True)
        
        assert isinstance(gbk_text, str)
        assert 'LOCUS' in gbk_text
        assert 'FEATURES' in gbk_text
        assert 'test_gene' in gbk_text
        assert sequence in gbk_text


class TestDataFrameProcessing:
    """Tests for DataFrame processing functions."""

    def test_get_clean_csv_df(self):
        """Test get_clean_csv_df function."""
        df = pd.DataFrame({
            'sseqid': ['gene1'],
            'qstart': [0],
            'qend': [10],
            'sframe': [1],
            'pident': [95.0],
            'slen': [100],
            'length': [10],
            'abs percmatch': [80.0],
            'fragment': [False],
            'db': ['test_db'],
            'Feature': ['Test Gene'],
            'Type': ['CDS'],
            'Description': ['Test description'],
            'qseq': ['ATCGATCGAT'],
            'extra_column': ['should_be_removed']
        })
        
        result = rsc.get_clean_csv_df(df)
        
        # Check that only expected columns are present
        expected_columns = [
            'sseqid', 'start location', 'end location', 'strand',
            'percent identity', 'full length of feature in db',
            'length of found feature', 'percent match length',
            'fragment', 'database', 'Feature', 'Type', 'Description', 'sequence'
        ]
        
        assert list(result.columns) == expected_columns
        assert 'extra_column' not in result.columns
        assert result['start location'].iloc[0] == 0
        assert result['database'].iloc[0] == 'test_db'


class TestConstants:
    """Tests for module constants."""

    def test_valid_extensions(self):
        """Test that valid extensions are defined correctly."""
        assert '.fasta' in rsc.valid_fasta_exts
        assert '.fa' in rsc.valid_fasta_exts
        assert '.fas' in rsc.valid_fasta_exts
        assert '.fna' in rsc.valid_fasta_exts
        
        assert '.gbk' in rsc.valid_genbank_exts
        assert '.gb' in rsc.valid_genbank_exts
        assert '.gbf' in rsc.valid_genbank_exts
        assert '.gbff' in rsc.valid_genbank_exts

    def test_max_plas_size(self):
        """Test that MAX_PLAS_SIZE is reasonable."""
        assert rsc.MAX_PLAS_SIZE == 50000
        assert isinstance(rsc.MAX_PLAS_SIZE, int)

    def test_df_cols(self):
        """Test that DF_COLS contains expected columns."""
        required_cols = [
            'sseqid', 'qstart', 'qend', 'sframe', 'score',
            'evalue', 'qseq', 'length', 'slen', 'pident'
        ]
        
        for col in required_cols:
            assert col in rsc.DF_COLS


class TestDatabaseHandling:
    """Tests for database caching utilities."""

    def test_get_db_dir_default(self, monkeypatch):
        monkeypatch.delenv("PLANNOTATE_DB_DIR", raising=False)
        monkeypatch.delenv("PLANNOTATE_SKIP_DB_DOWNLOAD", raising=False)
        expected = Path(user_cache_dir("pLannotate")) / "BLAST_dbs"
        assert rsc.get_db_dir(download=False) == expected

    def test_set_db_cache_dir(self, tmp_path, monkeypatch):
        monkeypatch.delenv("PLANNOTATE_SKIP_DB_DOWNLOAD", raising=False)
        rsc.set_db_cache_dir(tmp_path)
        assert rsc.get_db_dir(download=False) == tmp_path / "BLAST_dbs"

    def test_download_db(self, tmp_path, monkeypatch):
        monkeypatch.delenv("PLANNOTATE_SKIP_DB_DOWNLOAD", raising=False)
        monkeypatch.setenv("PLANNOTATE_AUTO_DOWNLOAD", "1")

        build_dir = tmp_path / "build"
        db_dir = build_dir / "BLAST_dbs"
        db_dir.mkdir(parents=True)
        (db_dir / "dummy.txt").write_text("hello")

        archive_path = tmp_path / "db.tar.gz"
        with tarfile.open(archive_path, "w:gz") as tar:
            tar.add(db_dir, arcname="BLAST_dbs")
        checksum = hashlib.sha256(archive_path.read_bytes()).hexdigest()

        result = rsc.download_db(cache_root=tmp_path, url=f"file://{archive_path}", checksum=checksum, force=True)
        assert (Path(result) / "dummy.txt").exists()


if __name__ == "__main__":
    pytest.main([__file__])