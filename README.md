# pLannotate-python

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python 3](https://img.shields.io/badge/Language-Python_3-steelblue.svg)
[![DOI](https://zenodo.org/badge/DOI/10.1093/nar/gkab374.svg)](https://doi.org/10.1093/nar/gkab374)

<img width="400" alt="pLannotate_logo" src="plannotate/data/images/pLannotate.png">

**Automated annotation of engineered plasmids**

`pLannotate-python` is a Python package for automatically annotating engineered plasmids using sequence similarity searches against curated databases. This is a streamlined, installable version of the original pLannotate tool, designed for programmatic use and integration into bioinformatics workflows. It has been optimized for performance with parallel processing.

## Features

- **Fast, parallel annotation**: Uses local installations of Diamond, BLAST, and Infernal to run comprehensive sequence searches concurrently.
- **Multiple databases**: Search against protein (fpbase, swissprot), nucleotide (snapgene), and RNA (Rfam) databases.
- **Circular plasmid support**: Handles origin-crossing features in circular plasmids.
- **Flexible output**: Generate GenBank files, CSV reports, or work with pandas DataFrames.
- **Batch processing**: Annotate multiple plasmids programmatically with high performance.

## Installation

### 1. Install pLannotate-python

```bash
# Install from PyPI
pip install plannotate-python

# Or install from source
git clone https://github.com/McClain-Thiel/pLannotate.git
cd pLannotate
pip install -e .

# Or install with uv (recommended)
uv pip install -e .
```

### 2. Install External Tools

`pLannotate-python` requires external bioinformatics tools for sequence searching.

**Required:**
- **BLAST+**: For nucleotide searches.
- **DIAMOND**: For fast protein searches.
- **Infernal**: For RNA secondary structure searches.
- **ripgrep (`rg`)**: For fast searches in compressed database files.

#### On macOS (using Homebrew)

```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install bioinformatics tools
brew install diamond blast infernal ripgrep
```

#### On Linux (using Conda/Mamba)

```bash
# Install conda/mamba if not already installed
# Then install bioinformatics tools
conda install -c bioconda diamond blast infernal ripgrep

# Or with mamba (faster)
mamba install -c bioconda diamond blast infernal ripgrep
```

#### On Linux (using package managers)

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install diamond-aligner ncbi-blast+ infernal ripgrep
```

**CentOS/RHEL/Fedora:**
```bash
# Install EPEL repository first
sudo yum install epel-release
sudo yum install diamond ncbi-blast+ infernal ripgrep
```

### 3. Verify Installation

```bash
# Check that tools are installed
diamond --version
blastn -version
cmscan -h
rg --version

# Test pLannotate import
python -c "from plannotate.annotate import annotate; print('✓ pLannotate-python installed successfully')"
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

## Database Setup

For full functionality, you need to set up local sequence databases.

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
# (You need to obtain the fpbase protein sequences in a FASTA file)
diamond makedb --in fpbase.fasta --db databases/fpbase

# Example: Create BLAST nucleotide database
makeblastdb -in snapgene.fasta -dbtype nucl -out databases/snapgene

# Example: Download and prepare Rfam (large download ~2GB)
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
mv Rfam.cm databases/
```

### 3. Update Database Configuration

Edit `plannotate/data/data/databases.yml` to point to your local database files. **Note:** You must provide the full path to the database files.

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

Rfam:
  method: infernal
  location: /path/to/your/databases/Rfam.cm
  priority: 3
  # ... other settings
```

## API Reference

### Core Functions

#### `annotate(sequence, yaml_file=None, linear=False, is_detailed=False)`

Annotate a DNA sequence by running searches against multiple databases in parallel.

**Parameters:**
- `sequence` (str): DNA sequence to annotate.
- `yaml_file` (str, optional): Path to a custom database configuration file. Defaults to the internal `databases.yml`.
- `linear` (bool): `True` for linear DNA, `False` for circular plasmids.
- `is_detailed` (bool): Include detailed feature information.

**Returns:**
- `pandas.DataFrame`: A DataFrame containing the annotation results.

#### `get_gbk(annotations_df, sequence, is_linear=False, record=None)`

Generate GenBank format output from an annotations DataFrame.

**Parameters:**
- `annotations_df` (DataFrame): Annotation results from `annotate()`.
- `sequence` (str): Original DNA sequence.
- `is_linear` (bool): `True` for linear DNA, `False` for circular.
- `record` (SeqRecord, optional): An existing Biopython SeqRecord to annotate.

**Returns:**
- `str`: GenBank formatted text.

## Troubleshooting

### Common Issues

**"Tool not found in PATH"**
```bash
# Ensure tools are installed and accessible in your shell's PATH
which diamond blastn cmscan rg

# If using conda, make sure your environment is activated
conda activate your_environment
```

**"No such file or directory" for databases**
- Verify that the database paths in `databases.yml` are correct and are absolute paths.
- Ensure the database files exist and have the correct read permissions.
- Check that Diamond databases have a `.dmnd` extension if you've specified one in the config.

**Empty results**
- The input sequence may not have any matches in the configured databases.
- Try lowering the identity thresholds or other search parameters in `databases.yml`.
- Verify that your databases are correctly formatted and contain relevant sequences for your plasmids.

### Performance
The annotation process is parallelized and will use multiple CPU cores to speed up the database searches. Performance will depend on the size of your databases and the number of available cores.

## Citation

If you use `pLannotate-python` in your research, please cite the original pLannotate paper:

> McGuffin, M.J., Thiel, M.C., Pineda, D.L. et al. pLannotate: automated annotation of engineered plasmids. *Nucleic Acids Research* (2021).

## License

This project is licensed under the GPL v3 License - see the `LICENSE` file for details.

## Links

- **Original pLannotate**: https://github.com/mmcguffi/pLannotate
- **Web server**: http://plannotate.barricklab.org/
- **This Fork**: https://github.com/McClain-Thiel/pLannotate
- **Issues**: https://github.com/McClain-Thiel/pLannotate/issues