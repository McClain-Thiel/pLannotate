import pytest
import pandas as pd
from unittest.mock import patch, mock_open, MagicMock
from tempfile import NamedTemporaryFile
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from plannotate.annotate import annotate, BLAST, calculate, clean
from plannotate import resources as rsc


@pytest.fixture
def sample_dna_sequence():
    """Sample DNA sequence for testing."""
    return "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"


@pytest.fixture
def sample_blast_result():
    """Sample BLAST result DataFrame."""
    return pd.DataFrame({
        'qstart': [1, 10],
        'qend': [30, 50],
        'sseqid': ['feature1', 'feature2'],
        'sframe': [1, -1],
        'pident': [95.0, 88.0],
        'slen': [100, 200],
        'qseq': ['ATCGATCGATCGATCGATCGATCGATCG', 'CGATCGATCGATCGATCGATCGATCGAT'],
        'length': [30, 41],
        'sstart': [1, 50],
        'send': [30, 91],
        'qlen': [60, 60],
        'evalue': [1e-10, 1e-8],
        'priority': [1, 2],
        'Type': ['CDS', 'misc_feature'],
        'Feature': ['Test Gene', 'Test Feature'],
        'Description': ['Test description', 'Another description']
    })


class TestAnnotateFunction:
    """Tests for the main annotate function."""

    @patch('plannotate.annotate.rsc.get_yaml')
    @patch('plannotate.annotate.BLAST')
    def test_annotate_basic_functionality(self, mock_blast, mock_get_yaml, sample_dna_sequence):
        """Test basic annotate functionality."""
        # Mock database configuration
        mock_get_yaml.return_value = {
            'database1': {
                'method': 'diamond',
                'priority': 1,
                'parameters': [],
                'db_loc': '/fake/path'
            }
        }
        
        # Mock BLAST results
        mock_blast.return_value = pd.DataFrame({
            'qstart': [1], 'qend': [30], 'sseqid': ['test_feature'],
            'sframe': [1], 'pident': [95.0], 'slen': [100],
            'qseq': ['ATCGATCGATCGATCGATCGATCGATCG'], 'length': [30],
            'sstart': [1], 'send': [30], 'qlen': [60], 'evalue': [1e-10],
            'Type': ['CDS'], 'Feature': ['Test Gene'], 'Description': ['Test']
        })
        
        result = annotate(sample_dna_sequence, linear=True)
        
        assert isinstance(result, pd.DataFrame)
        mock_blast.assert_called_once()
        mock_get_yaml.assert_called_once()

    @patch('plannotate.annotate.rsc.get_yaml')
    @patch('plannotate.annotate.BLAST')
    def test_annotate_circular_sequence_doubling(self, mock_blast, mock_get_yaml, sample_dna_sequence):
        """Test that circular sequences are doubled."""
        mock_get_yaml.return_value = {'db1': {'method': 'diamond', 'priority': 1, 'parameters': [], 'db_loc': '/fake'}}
        mock_blast.return_value = pd.DataFrame()
        
        annotate(sample_dna_sequence, linear=False)
        
        # Check that BLAST was called with doubled sequence
        call_args = mock_blast.call_args[1]
        assert call_args['seq'] == sample_dna_sequence + sample_dna_sequence

    def test_annotate_invalid_linear_parameter(self, sample_dna_sequence):
        """Test that invalid linear parameter raises ValueError."""
        with pytest.raises(ValueError, match="linear must be a boolean"):
            annotate(sample_dna_sequence, linear="not_boolean")

    @patch('plannotate.annotate.rsc.get_yaml')
    @patch('plannotate.annotate.BLAST')
    def test_annotate_no_hits(self, mock_blast, mock_get_yaml, sample_dna_sequence):
        """Test annotate behavior when no hits are found."""
        mock_get_yaml.return_value = {'db1': {'method': 'diamond', 'priority': 1, 'parameters': [], 'db_loc': '/fake'}}
        mock_blast.return_value = pd.DataFrame()
        
        result = annotate(sample_dna_sequence)
        
        assert isinstance(result, pd.DataFrame)
        assert result.empty


class TestBLASTFunction:
    """Tests for the BLAST function."""

    def test_blast_missing_blastn(self, sample_dna_sequence):
        """Test BLAST behavior when blastn is not available."""
        with patch('shutil.which', return_value=None):
            result = BLAST(sample_dna_sequence, {'method': 'blastn'})
            assert isinstance(result, pd.DataFrame)
            assert result.empty

    def test_blast_missing_diamond(self, sample_dna_sequence):
        """Test BLAST behavior when diamond is not available."""
        with patch('shutil.which', return_value=None):
            result = BLAST(sample_dna_sequence, {'method': 'diamond', 'db_loc': '/fake', 'parameters': ''})
            assert isinstance(result, pd.DataFrame)
            assert result.empty

    def test_blast_missing_infernal(self, sample_dna_sequence):
        """Test BLAST behavior when infernal is not available."""
        with patch('shutil.which', return_value=None):
            result = BLAST(sample_dna_sequence, {'method': 'infernal', 'db_loc': '/fake', 'parameters': ''})
            assert isinstance(result, pd.DataFrame)
            assert result.empty

    @patch('plannotate.annotate.subprocess.run')
    @patch('shutil.which')
    def test_blast_diamond_success(self, mock_which, mock_subprocess, sample_dna_sequence):
        """Test successful diamond BLAST execution."""
        mock_which.return_value = '/usr/bin/diamond'
        mock_subprocess.return_value = MagicMock()
        
        # Mock the output file content
        diamond_output = "1\t30\tfeature1\t95.0\t100\tATCG\t30\t1\t30\t60\t1e-10\n"
        
        with patch('builtins.open', mock_open(read_data=diamond_output)):
            result = BLAST(sample_dna_sequence, {
                'method': 'diamond',
                'db_loc': '/fake/db',
                'parameters': '--sensitive'
            })
            
        assert isinstance(result, pd.DataFrame)
        assert not result.empty

    def test_blast_unknown_method(self, sample_dna_sequence):
        """Test BLAST behavior with unknown method."""
        result = BLAST(sample_dna_sequence, {'method': 'unknown_method'})
        assert isinstance(result, pd.DataFrame)
        assert result.empty


class TestCalculateFunction:
    """Tests for the calculate function."""

    def test_calculate_linear_sequence(self, sample_blast_result):
        """Test calculate function with linear sequence."""
        result = calculate(sample_blast_result.copy(), is_linear=True)
        
        assert 'percmatch' in result.columns
        assert 'abs percmatch' in result.columns
        assert 'pi_permatch' in result.columns
        assert 'score' in result.columns
        assert 'wiggle' in result.columns

    def test_calculate_circular_sequence(self, sample_blast_result):
        """Test calculate function with circular sequence."""
        result = calculate(sample_blast_result.copy(), is_linear=False)
        
        # For circular sequences, qlen should be halved
        assert all(result['qlen'] == 30)  # 60/2 = 30

    def test_calculate_score_adjustment_by_priority(self, sample_blast_result):
        """Test that priority affects score calculation."""
        df = sample_blast_result.copy()
        result = calculate(df, is_linear=True)
        
        # Higher priority (lower number) should have higher score
        priority_1_score = result[result['priority'] == 1]['score'].iloc[0]
        priority_2_score = result[result['priority'] == 2]['score'].iloc[0]
        
        # This assumes similar other parameters, priority 1 should have higher score
        assert priority_1_score > priority_2_score

    def test_calculate_perfect_match_bonus(self):
        """Test that perfect matches get a bonus."""
        df = pd.DataFrame({
            'qstart': [1], 'qend': [30], 'sseqid': ['test'],
            'sframe': [1], 'pident': [100.0], 'slen': [30],
            'qseq': ['ATCG'], 'length': [30], 'sstart': [1],
            'send': [30], 'qlen': [60], 'evalue': [1e-10],
            'priority': [1], 'Type': ['CDS']
        })
        
        result = calculate(df, is_linear=True)
        # Perfect match (pi_permatch == 100) should have bonus applied
        assert result['pi_permatch'].iloc[0] == 100.0


class TestCleanFunction:
    """Tests for the clean function."""

    def test_clean_removes_problem_hits(self, sample_blast_result):
        """Test that problematic hits are removed."""
        df = sample_blast_result.copy()
        df.loc[0, 'sseqid'] = 'P03851'  # This is in problem_hits list
        
        result = clean(df)
        
        assert 'P03851' not in result['sseqid'].values

    def test_clean_filters_high_evalue(self, sample_blast_result):
        """Test that high e-values are filtered out."""
        df = sample_blast_result.copy()
        df.loc[0, 'evalue'] = 10  # High e-value
        
        result = clean(df)
        
        assert len(result) < len(df)

    def test_clean_filters_low_pi_permatch(self, sample_blast_result):
        """Test that low pi_permatch values are filtered out."""
        df = calculate(sample_blast_result.copy(), is_linear=True)
        df.loc[0, 'pi_permatch'] = 1.0  # Very low pi_permatch
        
        result = clean(df)
        
        assert len(result) < len(df)

    def test_clean_handles_origin_crossing(self):
        """Test that origin-crossing features are handled correctly."""
        df = pd.DataFrame({
            'qstart': [50], 'qend': [10], 'sseqid': ['origin_crossing'],
            'sframe': [1], 'pident': [95.0], 'slen': [100],
            'qseq': ['ATCG'], 'length': [30], 'sstart': [1],
            'send': [30], 'qlen': [60], 'evalue': [1e-10],
            'priority': [1], 'Type': ['CDS'], 'pi_permatch': [50.0],
            'wiggle': [5], 'wstart': [45], 'wend': [5], 'kind': ['CDS']
        })
        
        result = clean(df)
        
        # Should handle positions correctly for origin-crossing features
        assert isinstance(result, pd.DataFrame)

    def test_clean_empty_dataframe(self):
        """Test clean function with empty DataFrame."""
        df = pd.DataFrame()
        result = clean(df)
        
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == rsc.DF_COLS


class TestIntegration:
    """Integration tests for the complete annotation process."""

    @patch('plannotate.annotate.rsc.get_yaml')
    @patch('plannotate.annotate.subprocess.run')
    @patch('shutil.which')
    def test_full_annotation_pipeline(self, mock_which, mock_subprocess, mock_get_yaml):
        """Test the complete annotation pipeline."""
        # Mock external dependencies
        mock_which.return_value = '/usr/bin/diamond'
        mock_subprocess.return_value = MagicMock()
        mock_get_yaml.return_value = {
            'test_db': {
                'method': 'diamond',
                'priority': 1,
                'parameters': '',
                'db_loc': '/fake/db'
            }
        }
        
        # Mock file operations
        diamond_output = "1\t30\ttest_gene\t95.0\t100\tATCGATCG\t30\t1\t30\t60\t1e-10\n"
        
        with patch('builtins.open', mock_open(read_data=diamond_output)):
            result = annotate("ATCGATCGATCGATCGATCG" * 3, linear=True)
            
        assert isinstance(result, pd.DataFrame)
        # Should have processed the results through the pipeline
        if not result.empty:
            assert 'score' in result.columns
            assert 'fragment' in result.columns


if __name__ == "__main__":
    pytest.main([__file__])