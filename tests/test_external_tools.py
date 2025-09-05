import pytest
import pandas as pd
import os
import tempfile
import shutil
from pathlib import Path

from plannotate.annotate import BLAST
from plannotate import resources as rsc


class TestExternalTools:
    """Test pLannotate with actual external tools (diamond, blastn, cmscan)."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test databases."""
        with tempfile.TemporaryDirectory() as temp_dir:
            yield temp_dir

    @pytest.fixture
    def test_protein_sequence(self):
        """Test protein sequence for diamond testing."""
        # This is a short GFP-like sequence that might have matches
        return "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"

    @pytest.fixture
    def test_dna_sequence(self):
        """Test DNA sequence for blastn and cmscan testing."""
        # A sequence that includes some common features
        return ("ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGG"
                "GCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTATG"
                "GCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCA"
                "AGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCA"
                "ACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCC"
                "ACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCA"
                "CCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGT"
                "ACAAGTAA")

    def test_diamond_availability(self):
        """Test that diamond is available and working."""
        result = shutil.which("diamond")
        assert result is not None, "diamond not found in PATH"

    def test_blastn_availability(self):
        """Test that blastn is available and working."""
        result = shutil.which("blastn")
        assert result is not None, "blastn not found in PATH"

    def test_cmscan_availability(self):
        """Test that cmscan is available and working."""
        result = shutil.which("cmscan")
        assert result is not None, "cmscan not found in PATH"

    def test_create_diamond_database(self, temp_dir, test_protein_sequence):
        """Test creating a diamond database from a protein FASTA file."""
        # Create a test protein FASTA file
        protein_fasta = os.path.join(temp_dir, "test_proteins.fasta")
        with open(protein_fasta, 'w') as f:
            f.write(">test_protein_1|GFP\n")
            f.write("MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK\n")
            f.write(">test_protein_2|mCherry\n")
            f.write("MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK\n")

        # Create diamond database
        db_path = os.path.join(temp_dir, "test_proteins.dmnd")
        
        import subprocess
        result = subprocess.run([
            "diamond", "makedb", 
            "--in", protein_fasta,
            "--db", os.path.join(temp_dir, "test_proteins")
        ], capture_output=True, text=True)
        
        assert result.returncode == 0, f"Diamond makedb failed: {result.stderr}"
        assert os.path.exists(db_path), "Diamond database file not created"

    def test_create_blast_database(self, temp_dir, test_dna_sequence):
        """Test creating a BLAST database from a DNA FASTA file."""
        # Create a test DNA FASTA file
        dna_fasta = os.path.join(temp_dir, "test_dna.fasta")
        with open(dna_fasta, 'w') as f:
            f.write(">test_dna_1|GFP_gene\n")
            f.write(test_dna_sequence + "\n")
            f.write(">test_dna_2|promoter\n")
            f.write("TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC" * 10 + "\n")  # Synthetic promoter-like sequence

        # Create BLAST database
        import subprocess
        result = subprocess.run([
            "makeblastdb",
            "-in", dna_fasta,
            "-dbtype", "nucl",
            "-out", os.path.join(temp_dir, "test_dna")
        ], capture_output=True, text=True)
        
        assert result.returncode == 0, f"makeblastdb failed: {result.stderr}"

    def test_diamond_search_functionality(self, temp_dir, test_dna_sequence):
        """Test diamond functionality with BLAST function."""
        # Create minimal diamond database
        protein_fasta = os.path.join(temp_dir, "proteins.fasta")
        with open(protein_fasta, 'w') as f:
            f.write(">P42212|GFP_AEQVI Green fluorescent protein\n")
            f.write("MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK\n")

        import subprocess
        subprocess.run([
            "diamond", "makedb", 
            "--in", protein_fasta,
            "--db", os.path.join(temp_dir, "proteins")
        ], check=True)

        # Test diamond search
        db_config = {
            'method': 'diamond',
            'db_loc': os.path.join(temp_dir, "proteins.dmnd"),
            'parameters': '--id 50 --max-target-seqs 10'
        }
        
        result = BLAST(test_dna_sequence, db_config)
        
        assert isinstance(result, pd.DataFrame)
        # Diamond should find matches or return empty DataFrame (both are valid)

    def test_blastn_search_functionality(self, temp_dir, test_dna_sequence):
        """Test blastn functionality with BLAST function."""
        # Create minimal nucleotide database
        dna_fasta = os.path.join(temp_dir, "nucleotides.fasta")
        with open(dna_fasta, 'w') as f:
            f.write(">test_seq|GFP_gene\n")
            f.write(test_dna_sequence + "\n")

        import subprocess
        subprocess.run([
            "makeblastdb",
            "-in", dna_fasta,
            "-dbtype", "nucl",
            "-out", os.path.join(temp_dir, "nucleotides")
        ], check=True)

        # Test blastn search - this will use the online NCBI database
        db_config = {
            'method': 'blastn'
        }
        
        # Use a shorter sequence for online search to avoid timeout
        short_sequence = test_dna_sequence[:100]  # First 100 bp
        result = BLAST(short_sequence, db_config)
        
        assert isinstance(result, pd.DataFrame)
        # Online blastn might return results or empty DataFrame depending on sequence

    @pytest.mark.skip(reason="Requires Rfam database which is large to download")
    def test_cmscan_search_functionality(self, temp_dir, test_dna_sequence):
        """Test cmscan functionality with BLAST function."""
        # Note: This test is skipped because it requires the actual Rfam database
        # which is large (~2GB). In practice, users would need to download it.
        
        db_config = {
            'method': 'infernal',
            'db_loc': '/path/to/Rfam.cm',  # Would need actual Rfam database
            'parameters': '--cpu 1'
        }
        
        result = BLAST(test_dna_sequence, db_config)
        
        assert isinstance(result, pd.DataFrame)

    def test_blast_with_invalid_database(self, temp_dir):
        """Test BLAST functions handle invalid databases gracefully."""
        # Test diamond with non-existent database
        db_config = {
            'method': 'diamond',
            'db_loc': '/nonexistent/database.dmnd',
            'parameters': '--id 50'
        }
        
        result = BLAST("ATCGATCG", db_config)
        assert isinstance(result, pd.DataFrame)
        assert result.empty  # Should return empty DataFrame

    def test_tool_parameter_handling(self, temp_dir, test_dna_sequence):
        """Test that tool parameters are handled correctly."""
        # Create minimal diamond database
        protein_fasta = os.path.join(temp_dir, "test.fasta")
        with open(protein_fasta, 'w') as f:
            f.write(">test\nMVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLT\n")

        import subprocess
        subprocess.run([
            "diamond", "makedb", 
            "--in", protein_fasta,
            "--db", os.path.join(temp_dir, "test")
        ], check=True)

        # Test with various parameters
        db_config = {
            'method': 'diamond',
            'db_loc': os.path.join(temp_dir, "test.dmnd"),
            'parameters': '--id 90 --max-target-seqs 5 --sensitive'
        }
        
        result = BLAST(test_dna_sequence, db_config)
        assert isinstance(result, pd.DataFrame)


class TestRealWorldDatabases:
    """Test with databases that might exist on the system."""

    def test_existing_database_configuration(self):
        """Test that database configuration is loaded correctly."""
        yaml_path = rsc.get_yaml_path()
        assert os.path.exists(yaml_path)
        
        databases = rsc.get_yaml(yaml_path)
        assert isinstance(databases, dict)
        
        # Check expected databases are configured
        expected_dbs = ['fpbase', 'swissprot', 'snapgene', 'Rfam']
        for db in expected_dbs:
            assert db in databases
            assert 'method' in databases[db]
            assert 'priority' in databases[db]

    def test_database_metadata_files_exist(self):
        """Test that database metadata files exist."""
        data_dir = Path(__file__).parent.parent / "plannotate" / "data" / "data"
        
        # Check CSV files exist
        csv_files = ['fpbase.csv', 'snapgene.csv', 'swissprot.csv.gz']
        for csv_file in csv_files:
            csv_path = data_dir / csv_file
            assert csv_path.exists(), f"{csv_file} not found"

    def test_sample_sequences_exist(self):
        """Test that sample FASTA files exist for testing."""
        fastas_dir = Path(__file__).parent.parent / "plannotate" / "data" / "fastas"
        
        fasta_files = list(fastas_dir.glob("*.fa"))
        assert len(fasta_files) > 0, "No sample FASTA files found"
        
        # Test one sample file
        sample_fasta = fasta_files[0]
        with open(sample_fasta) as f:
            content = f.read()
            assert content.startswith(">"), "Invalid FASTA format"

    @pytest.mark.integration
    def test_annotation_with_sample_sequence(self):
        """Integration test with a sample sequence (requires external tool setup)."""
        from plannotate.annotate import annotate
        
        # Use a short test sequence
        test_seq = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGC"
        
        try:
            # This will use the default databases, which may not be fully set up
            # But should at least not crash
            result = annotate(test_seq, linear=True)
            assert isinstance(result, pd.DataFrame)
            # Result may be empty if databases aren't set up, but should not error
        except Exception as e:
            # If databases aren't configured, should get specific error messages
            assert "not found" in str(e) or "No such file" in str(e)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])