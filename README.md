# pLannotate-python

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python 3](https://img.shields.io/badge/Language-Python_3-steelblue.svg)
[![DOI](https://zenodo.org/badge/DOI/10.1093/nar/gkab374.svg)](https://doi.org/10.1093/nar/gkab374)

<img width="400" alt="pLannotate_logo" src="plannotate/data/images/pLannotate.png">

**Automated annotation of engineered plasmids**

`pLannotate-python` is a Python package for automatically annotating engineered plasmids using sequence similarity searches against curated databases. Fast, parallel processing with automatic database setup. All it is, is a python friendly wrapper around CLI tools. This means the CLI tools (and the databases they rely on) are required to be set up first. 

## Features

- **Fast, parallel annotation**: Uses Diamond, BLAST, and Infernal concurrently
- **Multiple databases**: Protein (fpbase, swissprot), nucleotide (snapgene), RNA (Rfam)
- **Circular plasmid support**: Handles origin-crossing features
- **Automatic database setup**: Downloads and configures databases (~900MB)
- **Flexible output**: GenBank files, CSV reports, or pandas DataFrames

## Installation

```bash
# Install with uv (recommended)
uv add plannotate-python

# Or with pip
pip install plannotate-python
```

### External Tools Required

```bash
# macOS (Homebrew)
brew install diamond blast infernal ripgrep

# Linux (conda/mamba)
conda install -c bioconda diamond blast infernal ripgrep

# Ubuntu/Debian
sudo apt install diamond-aligner ncbi-blast+ infernal ripgrep
```

### SSL Certificate Fix (macOS)
If you encounter SSL certificate errors during database download:
```bash
# Replace X.Y with your Python version (e.g., 3.11)
open "/Applications/Python X.Y/Install Certificates.command"
```

## "Quick" Start

**Automatic Database Setup:**
```python
import os
os.environ["PLANNOTATE_AUTO_DOWNLOAD"] = "1"  # Enable auto-download of databases
from plannotate.annotate import annotate

# First run will download databases (~900MB with progress bars)
>>> sequence="tgaccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttataggtctcaatccacgggtacgggtatggagaaacagtagagagttgcgataaaaagcgtcaggtagtatccgctaatcttatggataaaaatgctatggcatagcaaagtgtgacgccgtgcaaataatcaatgtggacttttctgccgtgattatagacacttttgttacgcgtttttgtcatggctttggtcccgctttgttacagaatgcttttaataagcggggttaccggtttggttagcgagaagagccagtaaaagacgcagtgacggcaatgtctgatgcaatatggacaattggtttcttgtaatcgttaatccgcaaataacgtaaaaacccgcttcggcgggtttttttatggggggagtttagggaaagagcatttgtcatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcgg"  # Your plasmid sequence
>>> result = annotate(sequence, linear=False)  # False for circular plasmids

>>> result
   qstart  qend              sseqid   pident  slen                                               qseq  length  ...  wiggle  wstart  wend  kind  qstart_dup qend_dup fragment
0     523   615   AmpR_promoter_(5)  100.000    92  TTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAG...      92  ...      13     536   601     1         523      614    False
1      11    83  rrnB_T1_terminator  100.000    72  CAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTA...      72  ...      10      21    72     1          11       82    False
2     155   440     araBAD_promoter   99.649   285  ATGGAGAAACAGTAGAGAGTTGCGATAAAAAGCGTCAGGTAGTATC...     285  ...      42     197   397     1         816     1100    False
3      98   126     T7Te_terminator  100.000    28                       GGCTCACCTTCGGGTGGGCCTTTCTGCG      28  ...       4     102   121     1          98      125    False
4     615   661                AmpR  100.000   861     ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGG      46  ...       6     621   654     1         615      660     True

[5 rows x 28 columns]


```

**Manual Database Setup:**
```python
from plannotate.resources import download_db
download_db()  # Downloads with progress bars and SSL error handling
```

**Generate GenBank Files:**
```python
from plannotate.resources import get_gbk
gbk_content = get_gbk(result, sequence, is_linear=False)
with open("my_plasmid.gbk", "w") as f:
    f.write(gbk_content)
```

## Configuration

**Environment Variables:**
- `PLANNOTATE_AUTO_DOWNLOAD=1` - Auto-download databases without prompting
- `PLANNOTATE_DB_DIR=/path` - Custom database directory
- `PLANNOTATE_SKIP_DB_DOWNLOAD=1` - Skip database downloads entirely

**Core Functions:**
- `annotate(sequence, linear=False)` - Annotate DNA sequence
- `get_gbk(annotations, sequence)` - Generate GenBank file
- `download_db()` - Download databases with progress bars

## Troubleshooting

**SSL Certificate Errors:** Run the SSL certificate fix command above
**Empty Results:** Sequence may not match database features  
**Tool Errors:** Ensure external tools are installed and in PATH

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