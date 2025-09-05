import pytest
import pandas as pd
from tempfile import NamedTemporaryFile
from unittest.mock import patch, mock_open
import os

from plannotate.infernal import parse_infernal


class TestParseInfernal:
    """Tests for the parse_infernal function."""

    def test_parse_infernal_basic(self):
        """Test basic infernal parsing functionality."""
        # Mock infernal output format
        infernal_output = """# cmscan :: search sequences against Rfam profiles
# INFERNAL 1.1.4 (Dec 2020)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#idx target name           accession      clan name            seq from    seq to        mdl from   mdl to       strand   trunc   pass   gc  bias  score   E-value inc description of target
#--- -------------------- -------------- -------------------- ----------- ----------- ---------- ---------- ---------- ------- ------ ---- ----- ------ --------- --- ---------------------
    1 5S_rRNA               RF00001        CL00113              100         200         1          100        +          no    1   0.45   0.0  100.5   1.2e-20 !   5S ribosomal RNA
    2 tRNA                  RF00005        CL00001              50          120         1          72         +          no    1   0.52   0.0   85.0   3.4e-15 !   transfer RNA
"""
        
        with NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp.write(infernal_output)
            tmp.flush()
            
            try:
                result = parse_infernal(tmp.name)
                
                assert isinstance(result, pd.DataFrame)
                assert len(result) == 2
                
                # Check column mapping
                assert 'sseqid' in result.columns
                assert 'qstart' in result.columns
                assert 'qend' in result.columns
                assert 'Feature' in result.columns
                assert 'Description' in result.columns
                
                # Check first row values
                first_row = result.iloc[0]
                assert first_row['sseqid'] == 1
                assert first_row['Feature'] == '5S rRNA'
                assert first_row['qstart'] == 99  # 100 - 1 (0-based)
                assert first_row['qend'] == 199   # 200 - 1 (0-based)
                assert first_row['sframe'] == 1
                assert first_row['pident'] == 100
                
            finally:
                os.unlink(tmp.name)

    def test_parse_infernal_empty_file(self):
        """Test parse_infernal with empty/no hits file."""
        infernal_output = """# cmscan :: search sequences against Rfam profiles
# INFERNAL 1.1.4 (Dec 2020)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#idx target name           accession      clan name            seq from    seq to        mdl from   mdl to       strand   trunc   pass   gc  bias  score   E-value inc description of target
#--- -------------------- -------------- -------------------- ----------- ----------- ---------- ---------- ---------- ------- ------ ---- ----- ------ --------- --- ---------------------
"""
        
        with NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp.write(infernal_output)
            tmp.flush()
            
            try:
                result = parse_infernal(tmp.name)
                
                assert isinstance(result, pd.DataFrame)
                assert result.empty
                
            finally:
                os.unlink(tmp.name)

    def test_parse_infernal_reverse_strand(self):
        """Test parse_infernal with reverse strand hit."""
        infernal_output = """# cmscan :: search sequences against Rfam profiles
# INFERNAL 1.1.4 (Dec 2020)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#idx target name           accession      clan name            seq from    seq to        mdl from   mdl to       strand   trunc   pass   gc  bias  score   E-value inc description of target
#--- -------------------- -------------- -------------------- ----------- ----------- ---------- ---------- ---------- ------- ------ ---- ----- ------ --------- --- ---------------------
    1 tRNA                  RF00005        CL00001              120         50          1          72         -          no    1   0.52   0.0   85.0   3.4e-15 !   transfer RNA
"""
        
        with NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp.write(infernal_output)
            tmp.flush()
            
            try:
                result = parse_infernal(tmp.name)
                
                assert len(result) == 1
                first_row = result.iloc[0]
                
                # Should swap start/end for reverse strand and adjust to 0-based
                assert first_row['qstart'] == 49   # min(120,50) - 1
                assert first_row['qend'] == 119    # max(120,50) - 1
                assert first_row['sframe'] == -1
                
            finally:
                os.unlink(tmp.name)

    def test_parse_infernal_column_processing(self):
        """Test that columns are processed correctly."""
        infernal_output = """# cmscan :: search sequences against Rfam profiles
# INFERNAL 1.1.4 (Dec 2020)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#idx target name           accession      clan name            seq from    seq to        mdl from   mdl to       strand   trunc   pass   gc  bias  score   E-value inc description of target
#--- -------------------- -------------- -------------------- ----------- ----------- ---------- ---------- ---------- ------- ------ ---- ----- ------ --------- --- ---------------------
    1 U6_snRNA              RF00026        CL00004              10          110         1          100        +          no    1   0.45   0.0  150.5   5.2e-25 !   U6 spliceosomal RNA
"""
        
        with NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp.write(infernal_output)
            tmp.flush()
            
            try:
                result = parse_infernal(tmp.name)
                
                assert len(result) == 1
                first_row = result.iloc[0]
                
                # Check feature name processing (underscores to spaces)
                assert first_row['Feature'] == 'U6 snRNA'
                
                # Check that accession and clan name dashes are replaced
                assert 'RF00026' in first_row['Description']
                assert 'CL00004' in first_row['Description']
                
                # Check calculated fields
                assert first_row['length'] == 101  # |qend - qstart| + 1
                assert first_row['slen'] == 100    # |send - sstart| + 1
                assert first_row['pident'] == 100  # Always set to 100 for infernal
                
                # Check that qseq is initialized as empty string
                assert first_row['qseq'] == ''
                
            finally:
                os.unlink(tmp.name)

    def test_parse_infernal_column_selection(self):
        """Test that only expected columns are returned."""
        infernal_output = """# cmscan :: search sequences against Rfam profiles
# INFERNAL 1.1.4 (Dec 2020)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#idx target name           accession      clan name            seq from    seq to        mdl from   mdl to       strand   trunc   pass   gc  bias  score   E-value inc description of target
#--- -------------------- -------------- -------------------- ----------- ----------- ---------- ---------- ---------- ------- ------ ---- ----- ------ --------- --- ---------------------
    1 5S_rRNA               RF00001        CL00113              100         200         1          100        +          no    1   0.45   0.0  100.5   1.2e-20 !   5S ribosomal RNA
"""
        
        with NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp.write(infernal_output)
            tmp.flush()
            
            try:
                result = parse_infernal(tmp.name)
                
                # Check that accession and clan name columns are dropped
                assert 'accession' not in result.columns
                assert 'clan name' not in result.columns
                
                # Check that expected columns are present
                expected_cols = [
                    'sseqid', 'qstart', 'qend', 'sstart', 'send',
                    'sframe', 'score', 'evalue', 'Feature', 'Description',
                    'qseq', 'length', 'slen', 'pident'
                ]
                
                for col in expected_cols:
                    assert col in result.columns
                    
            finally:
                os.unlink(tmp.name)

    def test_parse_infernal_numeric_conversion(self):
        """Test that numeric columns are properly converted."""
        infernal_output = """# cmscan :: search sequences against Rfam profiles
# INFERNAL 1.1.4 (Dec 2020)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#idx target name           accession      clan name            seq from    seq to        mdl from   mdl to       strand   trunc   pass   gc  bias  score   E-value inc description of target
#--- -------------------- -------------- -------------------- ----------- ----------- ---------- ---------- ---------- ------- ------ ---- ----- ------ --------- --- ---------------------
    1 5S_rRNA               RF00001        CL00113              100         200         5          95         +          no    1   0.45   0.0  100.5   1.2e-20 !   5S ribosomal RNA
"""
        
        with NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp.write(infernal_output)
            tmp.flush()
            
            try:
                result = parse_infernal(tmp.name)
                
                first_row = result.iloc[0]
                
                # Check that numeric columns have correct types
                assert isinstance(first_row['sseqid'], (int, pd.Int64Dtype))
                assert isinstance(first_row['qstart'], (int, pd.Int64Dtype))
                assert isinstance(first_row['qend'], (int, pd.Int64Dtype))
                assert isinstance(first_row['sstart'], (int, pd.Int64Dtype))
                assert isinstance(first_row['send'], (int, pd.Int64Dtype))
                assert isinstance(first_row['sframe'], (int, pd.Int64Dtype))
                assert isinstance(first_row['score'], (int, float, pd.Int64Dtype))
                assert isinstance(first_row['evalue'], (int, float, pd.Int64Dtype))
                assert isinstance(first_row['length'], (int, pd.Int64Dtype))
                assert isinstance(first_row['slen'], (int, pd.Int64Dtype))
                assert isinstance(first_row['pident'], (int, pd.Int64Dtype))
                
            finally:
                os.unlink(tmp.name)

    def test_parse_infernal_strand_replacement(self):
        """Test that strand symbols are correctly replaced."""
        infernal_output = """# cmscan :: search sequences against Rfam profiles
# INFERNAL 1.1.4 (Dec 2020)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#idx target name           accession      clan name            seq from    seq to        mdl from   mdl to       strand   trunc   pass   gc  bias  score   E-value inc description of target
#--- -------------------- -------------- -------------------- ----------- ----------- ---------- ---------- ---------- ------- ------ ---- ----- ------ --------- --- ---------------------
    1 5S_rRNA               RF00001        CL00113              100         200         1          100        +          no    1   0.45   0.0  100.5   1.2e-20 !   5S ribosomal RNA
    2 tRNA                  RF00005        CL00001              50          120         1          72         -          no    1   0.52   0.0   85.0   3.4e-15 !   transfer RNA
"""
        
        with NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp.write(infernal_output)
            tmp.flush()
            
            try:
                result = parse_infernal(tmp.name)
                
                # Check strand conversions: '+' -> 1, '-' -> -1
                assert result.iloc[0]['sframe'] == 1   # '+' strand
                assert result.iloc[1]['sframe'] == -1  # '-' strand
                
            finally:
                os.unlink(tmp.name)

    @patch('pandas.read_fwf')
    def test_parse_infernal_empty_data_error(self, mock_read_fwf):
        """Test parse_infernal handling of EmptyDataError."""
        mock_read_fwf.side_effect = pd.errors.EmptyDataError("No data")
        
        with NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp.write("# header only file\n")
            tmp.flush()
            
            try:
                result = parse_infernal(tmp.name)
                
                assert isinstance(result, pd.DataFrame)
                assert result.empty
                
            finally:
                os.unlink(tmp.name)

    def test_parse_infernal_malformed_file(self):
        """Test parse_infernal with malformed file."""
        # File without proper header structure
        malformed_output = "This is not a proper infernal output file"
        
        with NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp.write(malformed_output)
            tmp.flush()
            
            try:
                # Should not crash, but may return empty DataFrame
                result = parse_infernal(tmp.name)
                assert isinstance(result, pd.DataFrame)
                
            finally:
                os.unlink(tmp.name)


if __name__ == "__main__":
    pytest.main([__file__])