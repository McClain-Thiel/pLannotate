#!/usr/bin/env python3
"""
Test the full pLannotate workflow with external tools.
This script tests the complete annotation pipeline with a sample plasmid.
"""

import tempfile
import os
from pathlib import Path
from plannotate.annotate import annotate
from plannotate import resources as rsc


def create_test_diamond_db():
    """Create a minimal Diamond database for testing."""
    temp_dir = tempfile.mkdtemp()
    protein_fasta = os.path.join(temp_dir, "test_proteins.fasta")
    
    # Create a small protein database with common fluorescent proteins
    with open(protein_fasta, 'w') as f:
        f.write(">sp|P42212|GFP_AEQVI Enhanced green fluorescent protein\n")
        f.write("MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK\n")
        f.write(">sp|Q9U6Y8|MCHRY_DISCR mCherry fluorescent protein\n")
        f.write("MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK\n")
    
    # Build diamond database
    import subprocess
    db_path = os.path.join(temp_dir, "test_proteins")
    result = subprocess.run([
        "diamond", "makedb", 
        "--in", protein_fasta,
        "--db", db_path
    ], capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Diamond database creation failed: {result.stderr}")
        return None
        
    return db_path + ".dmnd"


def test_with_sample_sequence():
    """Test annotation with a sample GFP plasmid sequence."""
    
    print("=== Testing pLannotate with External Tools ===")
    print()
    
    # GFP coding sequence (shortened for testing)
    gfp_plasmid = (
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGG"
        "GCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTATG"
        "GCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCA"
        "AGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCA"
        "ACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCC"
        "ACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCA"
        "CCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGT"
        "ACAAGTAA"
    )
    
    print(f"Test sequence length: {len(gfp_plasmid)} bp")
    print()
    
    # Test 1: With tool availability warnings (no databases)
    print("1. Testing with default configuration (expected warnings about missing databases):")
    try:
        result = annotate(gfp_plasmid, linear=True)
        print(f"   Result: {type(result)}, Shape: {result.shape}, Empty: {result.empty}")
        if not result.empty:
            print(f"   Found {len(result)} annotations")
            print("   Features:", result['Feature'].tolist() if 'Feature' in result.columns else 'N/A')
    except Exception as e:
        print(f"   Expected result - no databases configured: {e}")
    print()
    
    # Test 2: With a minimal diamond database
    print("2. Testing with minimal Diamond database:")
    diamond_db = create_test_diamond_db()
    
    if diamond_db:
        # Create temporary YAML config
        import yaml
        import tempfile
        
        temp_yaml = tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False)
        test_config = {
            'test_fpbase': {
                'version': 'test',
                'method': 'diamond',
                'location': diamond_db,
                'priority': 1,
                'parameters': ['--id', '50', '--max-target-seqs', '10'],
                'details': {
                    'default_type': 'CDS',
                    'location': 'None',
                    'compressed': False
                }
            }
        }
        
        yaml.dump(test_config, temp_yaml)
        temp_yaml.close()
        
        try:
            result = annotate(gfp_plasmid, yaml_file=temp_yaml.name, linear=True)
            print(f"   Result: {type(result)}, Shape: {result.shape}, Empty: {result.empty}")
            
            if not result.empty:
                print(f"   Found {len(result)} annotations")
                for i, row in result.iterrows():
                    feature = row.get('Feature', 'Unknown')
                    start = row.get('qstart', 'N/A')
                    end = row.get('qend', 'N/A') 
                    score = row.get('score', 'N/A')
                    print(f"   - {feature}: {start}-{end} (score: {score})")
            else:
                print("   No annotations found (expected - sequence may be too short for matches)")
                
        except Exception as e:
            print(f"   Error: {e}")
        finally:
            os.unlink(temp_yaml.name)
            # Clean up database files
            db_dir = os.path.dirname(diamond_db)
            import shutil
            shutil.rmtree(db_dir)
    else:
        print("   Could not create Diamond database")
    print()
    
    # Test 3: Tool availability
    print("3. External tool availability check:")
    import shutil
    tools = ['diamond', 'blastn', 'cmscan']
    for tool in tools:
        available = shutil.which(tool) is not None
        print(f"   {tool}: {'✓ Available' if available else '✗ Not found'}")
    print()
    
    # Test 4: Sample FASTA files
    print("4. Sample FASTA files in repo:")
    fastas_dir = Path(__file__).parent / "plannotate" / "data" / "fastas"
    if fastas_dir.exists():
        fasta_files = list(fastas_dir.glob("*.fa"))[:5]  # Show first 5
        for fasta in fasta_files:
            print(f"   - {fasta.name}")
        if len(fasta_files) >= 5:
            print("   - ... and more")
    else:
        print("   No FASTA directory found")
    print()
    
    print("=== Test Complete ===")


if __name__ == "__main__":
    test_with_sample_sequence()