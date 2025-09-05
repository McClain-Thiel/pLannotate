# pLannotate

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python 3](https://img.shields.io/badge/Language-Python_3-steelblue.svg)
[![DOI](https://zenodo.org/badge/DOI/10.1093/nar/gkab374.svg)](https://doi.org/10.1093/nar/gkab374)

<img width="400" alt="pLannotate_logo" src="plannotate/data/images/pLannotate.png">

**Automated annotation of engineered plasmids**

pLannotate is a Python package for automatically annotating engineered plasmids using sequence similarity searches against curated databases. This is a streamlined, installable version of the original pLannotate tool, designed for programmatic use and integration into bioinformatics workflows.

## Features

- **Fast annotation**: Uses Diamond, BLAST, and Infernal for comprehensive sequence searches
- **Multiple databases**: Search against protein (fpbase, swissprot), nucleotide (snapgene), and RNA (Rfam) databases
- **Circular plasmid support**: Handles origin-crossing features in circular plasmids
- **Flexible output**: Generate GenBank files, CSV reports, or work with pandas DataFrames
- **Batch processing**: Annotate multiple plasmids programmatically

## Installation

### 1. Install pLannotate

```bash
# Install from PyPI (when available)
pip install plannotate

# Or install from source
git clone https://github.com/McClain-Thiel/pLannotate.git
cd pLannotate
pip install -e .

# Or install with uv (recommended)
uv add .
```

### 2. Install External Tools

pLannotate requires external bioinformatics tools for sequence searching:

#### On macOS (using Homebrew)

```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install bioinformatics tools
brew install diamond
brew install blast
brew install infernal
```

#### On Linux (using Conda/Mamba)

```bash
# Install conda/mamba if not already installed
# Then install bioinformatics tools
conda install -c bioconda diamond blast infernal

# Or with mamba (faster)
mamba install -c bioconda diamond blast infernal
```

#### On Linux (using package managers)

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install diamond-aligner ncbi-blast+ infernal
```

**CentOS/RHEL/Fedora:**
```bash
# Install EPEL repository first
sudo yum install epel-release
sudo yum install diamond ncbi-blast+ infernal
```

### 3. Verify Installation

```bash
# Check that tools are installed
diamond version
blastn -version
cmscan -h

# Test pLannotate import
python -c "from plannotate.annotate import annotate; print('✓ pLannotate installed successfully')"
```

## Quick Start

### Basic Usage

```python
from plannotate.annotate import annotate

# Annotate a plasmid sequence
sequence = "ATGGTGAGCAAGGGCGAGGAGCTG..."  # Your plasmid sequence
result = annotate(sequence, linear=False)  # False for circular plasmids

# View results
print(f"Found {len(result)} annotations")
print(result[['Feature', 'Type', 'qstart', 'qend', 'pident']].head())
```

### Generate GenBank File

```python
from plannotate.annotate import annotate
from plannotate.resources import get_gbk

# Annotate sequence
sequence = "ATGGTGAGCAAGGGCGAGGAGCTG..."
annotations = annotate(sequence, linear=False)

# Generate GenBank file
gbk_content = get_gbk(annotations, sequence, is_linear=False)

# Save to file
with open("my_plasmid.gbk", "w") as f:
    f.write(gbk_content)
```

### Working with Sample Plasmids

```python
from pathlib import Path
from Bio import SeqIO
from plannotate.annotate import annotate

# Use included sample plasmids
sample_dir = Path("plannotate/data/fastas")
for fasta_file in sample_dir.glob("*.fa"):
    # Load sequence
    record = list(SeqIO.parse(fasta_file, "fasta"))[0]
    sequence = str(record.seq)
    
    # Annotate
    result = annotate(sequence, linear=False)
    print(f"{fasta_file.name}: {len(result)} features found")
```

## Database Setup

For full functionality, you need to set up sequence databases:

### 1. Download/Create Databases

**Protein Databases (Diamond format):**
- **fpbase**: Fluorescent proteins database
- **swissprot**: SwissProt protein database

**Nucleotide Databases (BLAST format):**
- **snapgene**: Common cloning features

**RNA Databases (Infernal format):**
- **Rfam**: RNA families database

### 2. Example Database Setup

```bash
# Create database directory
mkdir -p databases

# Example: Create fpbase diamond database
# (You need to obtain the fpbase protein sequences)
diamond makedb --in fpbase.fasta --db databases/fpbase

# Example: Create BLAST nucleotide database
makeblastdb -in snapgene.fasta -dbtype nucl -out databases/snapgene

# Example: Download and prepare Rfam (large download ~2GB)
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
mv Rfam.cm databases/
```

### 3. Update Database Configuration

Edit `plannotate/data/data/databases.yml` to point to your database files:

```yaml
fpbase:
  method: diamond
  location: /path/to/your/databases/fpbase.dmnd
  priority: 1
  # ... other settings

snapgene:
  method: blastn
  location: /path/to/your/databases/snapgene
  priority: 1
  # ... other settings
```

## Advanced Usage

### Custom Database Configuration

```python
# Use custom database configuration
custom_config = "my_databases.yml"
result = annotate(sequence, yaml_file=custom_config, linear=False)
```

### Batch Processing

```python
import pandas as pd
from plannotate.annotate import annotate

sequences = {
    "plasmid1": "ATGGTGAGCAAG...",
    "plasmid2": "ATGGTGAGCAAG...",
    # ... more sequences
}

results = []
for name, seq in sequences.items():
    annotations = annotate(seq, linear=False)
    annotations['plasmid_name'] = name
    results.append(annotations)

# Combine all results
all_annotations = pd.concat(results, ignore_index=True)
all_annotations.to_csv("batch_annotations.csv", index=False)
```

### Filter Results

```python
# Get only CDS features with high identity
cds_features = result[
    (result['Type'] == 'CDS') & 
    (result['pident'] > 90)
]

# Get features above a certain score threshold
high_score_features = result[result['score'] > 100]
```

## API Reference

### Core Functions

#### `annotate(sequence, yaml_file=None, linear=False, is_detailed=False)`

Annotate a DNA sequence.

**Parameters:**
- `sequence` (str): DNA sequence to annotate
- `yaml_file` (str, optional): Path to database configuration file
- `linear` (bool): True for linear DNA, False for circular plasmids
- `is_detailed` (bool): Include detailed feature information

**Returns:**
- `pandas.DataFrame`: Annotation results

#### `get_gbk(annotations_df, sequence, is_linear=False, record=None)`

Generate GenBank format output.

**Parameters:**
- `annotations_df` (DataFrame): Annotation results from `annotate()`
- `sequence` (str): Original DNA sequence
- `is_linear` (bool): True for linear DNA, False for circular
- `record` (SeqRecord, optional): Existing SeqRecord to annotate

**Returns:**
- `str`: GenBank formatted text

### DataFrame Columns

The annotation results DataFrame contains these key columns:

- `Feature`: Feature name/description
- `Type`: Feature type (CDS, misc_feature, etc.)
- `qstart`, `qend`: Start and end positions (0-based)
- `pident`: Percent identity
- `length`: Feature length
- `score`: Annotation confidence score
- `fragment`: Boolean indicating if feature is truncated
- `db`: Source database

## Troubleshooting

### Common Issues

**"Tool not found in PATH"**
```bash
# Ensure tools are installed and accessible
which diamond blastn cmscan

# If using conda, activate the environment
conda activate your_environment
```

**"No such file or directory" for databases**
- Verify database paths in `databases.yml`
- Ensure database files exist and have correct permissions
- Check that Diamond databases have `.dmnd` extension

**Empty results**
- Sequence may not have matches in current databases
- Try lowering identity thresholds in database parameters
- Verify databases contain relevant sequences for your plasmids

### Performance Tips

- Use smaller, curated databases for faster searches
- Adjust database parameters (identity thresholds, max targets)
- For batch processing, consider parallel execution

## Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## Citation

If you use pLannotate in your research, please cite:

> Barrick Lab. pLannotate: automated annotation of engineered plasmids. *Nucleic Acids Research* (2021).

## License

This project is licensed under the GPL v3 License - see the LICENSE file for details.

## Links

- Original pLannotate: https://github.com/mmcguffi/pLannotate
- Web server: http://plannotate.barricklab.org/
- Documentation: https://github.com/McClain-Thiel/pLannotate#readme
- Issues: https://github.com/McClain-Thiel/pLannotate/issues