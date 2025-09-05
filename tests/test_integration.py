import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from tempfile import NamedTemporaryFile
import os

from plannotate.annotate import annotate, BLAST, calculate, clean
from plannotate import resources as rsc


class TestBasicFunctionality:
    """Integration tests to verify the package works as expected."""

    def test_package_import(self):
        """Test that the package can be imported successfully."""
        from plannotate.annotate import annotate
        from plannotate import resources
        from plannotate.infernal import parse_infernal
        from plannotate.blast_parser import _parse_blast_xml
        assert callable(annotate)

    @patch('plannotate.annotate.rsc.get_yaml')
    @patch('plannotate.annotate.BLAST')
    def test_annotate_with_mocked_dependencies(self, mock_blast, mock_get_yaml):
        """Test annotate function with mocked external dependencies."""
        # Setup mocks
        mock_get_yaml.return_value = {
            'test_db': {
                'method': 'diamond',
                'priority': 1,
                'parameters': ['--sensitive'],
                'db_loc': '/fake/db'
            }
        }
        
        mock_blast_result = pd.DataFrame({
            'qstart': [10, 50],
            'qend': [40, 100],
            'sseqid': ['gene1', 'gene2'],
            'sframe': [1, -1],
            'pident': [95.0, 88.0],
            'slen': [100, 200],
            'qseq': ['ATCGATCGATCGATCGATCGATCGATCG', 'CGATCGATCGATCGATCGATCGATCGAT'],
            'length': [31, 51],
            'sstart': [1, 50],
            'send': [31, 101],
            'qlen': [200, 200],
            'evalue': [1e-20, 1e-15],
            'Type': ['CDS', 'misc_feature'],
            'Feature': ['Test Gene 1', 'Test Feature 2'],
            'Description': ['Description 1', 'Description 2']
        })
        mock_blast.return_value = mock_blast_result
        
        # Test the annotate function
        sequence = "ATCGATCGATCGATCGATCG" * 10  # 200 bp sequence
        result = annotate(sequence, linear=True)
        
        # Verify results
        assert isinstance(result, pd.DataFrame)
        assert not result.empty
        assert len(result) <= 2  # May be filtered during processing
        
        # Check that processing occurred
        if not result.empty:
            assert 'score' in result.columns
            assert 'fragment' in result.columns
            assert 'percmatch' in result.columns

    @patch('plannotate.annotate.rsc.get_yaml')
    @patch('plannotate.annotate.BLAST')
    def test_annotate_circular_vs_linear(self, mock_blast, mock_get_yaml):
        """Test that circular and linear modes produce different query sequences."""
        mock_get_yaml.return_value = {
            'test_db': {
                'method': 'diamond',
                'priority': 1,
                'parameters': [],
                'db_loc': '/fake/db'
            }
        }
        
        mock_blast.return_value = pd.DataFrame()  # No hits
        
        sequence = "ATCGATCGATCGATCG"
        
        # Test linear
        annotate(sequence, linear=True)
        linear_call_args = mock_blast.call_args
        
        # Reset mock
        mock_blast.reset_mock()
        
        # Test circular
        annotate(sequence, linear=False)
        circular_call_args = mock_blast.call_args
        
        # For circular, sequence should be doubled
        assert len(circular_call_args[1]['seq']) == 2 * len(linear_call_args[1]['seq'])

    def test_calculate_function_basic(self):
        """Test the calculate function works correctly."""
        test_df = pd.DataFrame({
            'qstart': [10, 50],
            'qend': [40, 100],
            'sseqid': ['gene1', 'gene2'],
            'sframe': [1, -1],
            'pident': [95.0, 88.0],
            'slen': [100, 200],
            'qseq': ['ATCG' * 7 + 'ATC', 'CGAT' * 12 + 'CGT'],
            'length': [31, 51],
            'sstart': [1, 50],
            'send': [31, 101],
            'qlen': [200, 200],
            'evalue': [1e-20, 1e-15],
            'priority': [1, 2],
            'Type': ['CDS', 'misc_feature']
        })
        
        result = calculate(test_df, is_linear=True)
        
        # Check that new columns are added
        expected_cols = ['percmatch', 'abs percmatch', 'pi_permatch', 'score', 'wiggle', 'wstart', 'wend']
        for col in expected_cols:
            assert col in result.columns
        
        # Check that positions are adjusted (0-based)
        assert result.iloc[0]['qstart'] == 9  # 10 - 1
        assert result.iloc[0]['qend'] == 39   # 40 - 1

    def test_clean_function_basic(self):
        """Test the clean function works correctly."""
        test_df = pd.DataFrame({
            'qstart': [10, 50],
            'qend': [40, 100],
            'sseqid': ['gene1', 'gene2'],
            'sframe': [1, -1],
            'pident': [95.0, 88.0],
            'slen': [100, 200],
            'qseq': ['ATCG' * 7 + 'ATC', 'CGAT' * 12 + 'CGT'],
            'length': [31, 51],
            'sstart': [1, 50],
            'send': [31, 101],
            'qlen': [200, 200],
            'evalue': [1e-20, 1e-15],
            'priority': [1, 2],
            'Type': ['CDS', 'misc_feature'],
            'percmatch': [31.0, 25.5],
            'abs percmatch': [31.0, 25.5],
            'pi_permatch': [29.45, 22.44],  # Above 3% threshold
            'wiggle': [4, 7],
            'wstart': [14, 57],
            'wend': [36, 93],
            'kind': [1, 1]
        })
        
        result = clean(test_df)
        
        # Should return a DataFrame
        assert isinstance(result, pd.DataFrame)
        
        # Check that duplicates are removed and DataFrame is reset
        assert result.index.tolist() == list(range(len(result)))

    def test_resources_basic_functions(self):
        """Test basic resource functions work."""
        # Test path functions don't crash
        yaml_path = rsc.get_yaml_path()
        assert isinstance(yaml_path, str)
        assert 'databases.yml' in yaml_path
        
        name, ext = rsc.get_name_ext('test.fasta')
        assert name == 'test'
        assert ext == '.fasta'
        
        # Test sequence validation
        rsc.validate_sequence('ATCGATCG')  # Should not raise
        
        with pytest.raises(ValueError):
            rsc.validate_sequence('ATCGXYZ')  # Invalid characters

    @patch('plannotate.annotate.subprocess.run')
    @patch('shutil.which')
    def test_blast_tool_availability_handling(self, mock_which, mock_subprocess):
        """Test that BLAST handles missing external tools gracefully."""
        # Test when tools are not available
        mock_which.return_value = None
        
        result_diamond = BLAST("ATCG", {'method': 'diamond', 'db_loc': '/fake', 'parameters': ''})
        result_infernal = BLAST("ATCG", {'method': 'infernal', 'db_loc': '/fake', 'parameters': ''})
        result_blastn = BLAST("ATCG", {'method': 'blastn'})
        
        # Should return empty DataFrames with correct columns
        assert isinstance(result_diamond, pd.DataFrame)
        assert isinstance(result_infernal, pd.DataFrame)
        assert isinstance(result_blastn, pd.DataFrame)
        
        assert result_diamond.empty
        assert result_infernal.empty
        assert result_blastn.empty

    def test_yaml_parsing(self):
        """Test YAML configuration parsing."""
        yaml_content = """
test_db:
  method: diamond
  location: /path/to/db
  priority: 1
  parameters:
    - "--sensitive"
    - "--max-target-seqs"
    - "500"
"""
        with NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as tmp:
            tmp.write(yaml_content)
            tmp.flush()
            
            try:
                result = rsc.get_yaml(tmp.name)
                assert 'test_db' in result
                assert result['test_db']['method'] == 'diamond'
                assert result['test_db']['db_loc'] == '/path/to/db'
                assert result['test_db']['parameters'] == '--sensitive --max-target-seqs 500'
            finally:
                os.unlink(tmp.name)

    def test_genbank_generation_basic(self):
        """Test basic GenBank file generation."""
        test_df = pd.DataFrame({
            'qstart': [0],
            'qend': [10],
            'qlen': [20],
            'sframe': [1],
            'Feature': ['test_gene'],
            'Type': ['CDS'],
            'db': ['test_db'],
            'pident': [95.0],
            'percmatch': [50.0],
            'fragment': [False]
        })
        
        sequence = "ATCGATCGATCGATCGATCG"
        
        # Test SeqRecord generation
        record = rsc.get_seq_record(test_df, sequence, is_linear=True)
        assert str(record.seq) == sequence
        assert len(record.features) == 1
        assert record.features[0].qualifiers['label'][0] == 'test_gene'
        
        # Test GenBank text generation
        gbk_text = rsc.get_gbk(test_df, sequence, is_linear=True)
        assert isinstance(gbk_text, str)
        assert 'LOCUS' in gbk_text
        assert 'test_gene' in gbk_text

    @patch('plannotate.annotate.rsc.get_yaml')
    def test_annotate_error_handling(self, mock_get_yaml):
        """Test error handling in annotate function."""
        mock_get_yaml.return_value = {}
        
        # Test with invalid linear parameter
        with pytest.raises(ValueError, match="linear must be a boolean"):
            annotate("ATCG", linear="invalid")
        
        # Test with empty databases (should return empty DataFrame)
        result = annotate("ATCG", linear=True)
        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_module_constants(self):
        """Test that module constants are properly defined."""
        # Test file extensions
        assert '.fasta' in rsc.valid_fasta_exts
        assert '.gbk' in rsc.valid_genbank_exts
        
        # Test max plasmid size
        assert isinstance(rsc.MAX_PLAS_SIZE, int)
        assert rsc.MAX_PLAS_SIZE > 0
        
        # Test DataFrame columns
        assert isinstance(rsc.DF_COLS, list)
        assert 'sseqid' in rsc.DF_COLS
        assert 'qstart' in rsc.DF_COLS


class TestRealWorldScenarios:
    """Test realistic usage scenarios."""

    @patch('plannotate.annotate.rsc.get_yaml')
    @patch('plannotate.annotate.BLAST')
    def test_typical_plasmid_annotation_workflow(self, mock_blast, mock_get_yaml):
        """Test a typical plasmid annotation workflow."""
        # Mock databases similar to real configuration
        mock_get_yaml.return_value = {
            'fpbase': {
                'method': 'diamond',
                'priority': 1,
                'parameters': ['-k', '0', '--sensitive'],
                'db_loc': '/fake/fpbase'
            },
            'rfam': {
                'method': 'infernal', 
                'priority': 3,
                'parameters': ['--cpu', '1'],
                'db_loc': '/fake/rfam'
            }
        }
        
        # Mock realistic BLAST results
        mock_blast.side_effect = [
            # FPbase results (proteins)
            pd.DataFrame({
                'qstart': [100, 500],
                'qend': [800, 1200],
                'sseqid': ['mCherry', 'EGFP'],
                'sframe': [1, 1],
                'pident': [98.5, 95.2],
                'slen': [700, 700],
                'qseq': ['ATG' + 'N' * 697, 'ATG' + 'N' * 697],
                'length': [701, 701],
                'sstart': [1, 1],
                'send': [700, 700],
                'qlen': [3000, 3000],
                'evalue': [1e-100, 1e-95],
                'Type': ['CDS', 'CDS'],
                'Feature': ['mCherry fluorescent protein', 'Enhanced GFP'],
                'Description': ['Red fluorescent protein', 'Green fluorescent protein']
            }),
            # Rfam results (RNAs)
            pd.DataFrame({
                'qstart': [2000],
                'qend': [2080],
                'sseqid': ['tRNA-Ala'],
                'sframe': [1],
                'pident': [100.0],
                'slen': [80],
                'qseq': ['N' * 80],
                'length': [80],
                'sstart': [1],
                'send': [80],
                'qlen': [3000],
                'evalue': [1e-40],
                'Type': ['tRNA'],
                'Feature': ['tRNA-Ala'],
                'Description': ['Alanine transfer RNA']
            })
        ]
        
        # Test with a typical plasmid sequence
        plasmid_seq = "ATCG" * 750  # 3000 bp plasmid
        
        result = annotate(plasmid_seq, linear=False)  # Circular plasmid
        
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            # Should have processed features with scores, fragments, etc.
            expected_cols = ['Feature', 'Type', 'qstart', 'qend', 'score', 'fragment']
            for col in expected_cols:
                assert col in result.columns
            
            # Features should be properly categorized
            if 'Type' in result.columns:
                types = result['Type'].unique()
                assert all(t in ['CDS', 'tRNA', 'misc_feature'] for t in types)

    def test_edge_cases(self):
        """Test edge cases and boundary conditions."""
        # Very short sequence
        short_seq = "ATCG"
        rsc.validate_sequence(short_seq)  # Should not raise
        
        # Sequence at max length
        max_seq = "A" * rsc.MAX_PLAS_SIZE
        rsc.validate_sequence(max_seq)  # Should not raise
        
        # Sequence just over max length
        with pytest.raises(ValueError):
            rsc.validate_sequence("A" * (rsc.MAX_PLAS_SIZE + 1))


if __name__ == "__main__":
    pytest.main([__file__])