# ROSE2: Rank Ordering of Super-Enhancers

[![PyPI version](https://badge.fury.io/py/rose2.svg)](https://badge.fury.io/py/rose2)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

ROSE (Rank Ordering of Super-Enhancers) is a computational pipeline for the identification and analysis of super-enhancers from ChIP-seq data.

## Overview

ROSE enables researchers to:
- Identify super-enhancers from ChIP-seq data
- Rank enhancer regions by signal density
- Map stitched enhancers to genes
- Analyze regulatory landscapes

## Installation

### From PyPI (Recommended)

```bash
pip install rose2
```

### From Source

```bash
git clone https://github.com/stjude/ROSE2.git
cd ROSE2
pip install -e .
```

## Dependencies

- Python 3.8 or higher
- [samtools](http://www.htslib.org/) (for BAM file processing)
- R (version ≥ 3.4) (for statistical analysis and plotting)
- [bedtools](https://bedtools.readthedocs.io/) (version ≥ 2)

Install system dependencies on Ubuntu/Debian:
```bash
sudo apt-get install samtools bedtools r-base
```

On macOS with Homebrew:
```bash
brew install samtools bedtools r
```

## Quick Start

### Full ROSE Analysis

```bash
rose2 main \
    -i input_peaks.gff \
    -r sample.bam \
    -g HG19 \
    -o output_directory/
```

### Map BAM to GFF Regions

```bash
rose2 bamToGFF \
    -b sample.bam \
    -i regions.gff \
    -m 50 \
    -o output.gff
```

### Map Enhancers to Genes

```bash
rose2 geneMapper \
    -i enhancers_table.txt \
    -g HG19 \
    -o output_directory/
```

## Usage

### Main ROSE Pipeline

The main command runs the full ROSE analysis:

```bash
rose2 main -g [GENOME] -i [INPUT_GFF] -r [RANKBY_BAM] -o [OUTPUT_FOLDER] [OPTIONS]
```

**Required Arguments:**
- `-i, --input`: Input .gff or .bed file of binding sites
- `-r, --rankby`: BAM file to rank enhancers by
- `-o, --out`: Output folder
- `-g, --genome`: Genome build (HG18, HG19, HG38, MM8, MM9, MM10) *OR*
- `--custom`: Path to custom genome annotation file

**Optional Arguments:**
- `-b, --bams`: Comma-separated list of additional BAM files
- `-c, --control`: Control BAM file
- `-s, --stitch`: Maximum linking distance for stitching (default: 12500)
- `-t, --tss`: Distance from TSS to exclude (default: 0, no exclusion)

### BAM to GFF Mapping

Map read density from BAM files to genomic regions:

```bash
rose2 bamToGFF -b [BAM_FILE] -i [GFF_FILE] -o [OUTPUT] [OPTIONS]
```

**Options:**
- `-s, --sense`: Map to '+', '-', or 'both' strands (default: both)
- `-e, --extension`: Extend reads by N bp (default: 200)
- `-r, --rpm`: Normalize to reads per million
- `-m, --matrix`: Number of bins for matrix output

### Gene Mapper

Map enhancers to nearby genes:

```bash
rose2 geneMapper -g [GENOME] -i [ENHANCER_FILE] [OPTIONS]
```

**Options:**
- `-l, --list`: Gene list to filter
- `-o, --out`: Output folder
- `-r, --refseq`: Output by RefSeq ID instead of gene name
- `-c, --control`: Subtract input signal from sample

## Input Requirements

### BAM Files
- Must be coordinate-sorted
- Must be indexed (`.bai` file present)
- Chromosome names must start with `chr`

To sort and index:
```bash
samtools sort input.bam -o sorted.bam
samtools index sorted.bam
```

### Peak/Region Files (GFF or BED format)

**GFF Format:**
```
chr1    peak_1    .    1000    2000    .    +    .    peak_1
chr1    peak_2    .    3000    4000    .    +    .    peak_2
```

**BED Format:**
```
chr1    1000    2000    peak_1    100    +
chr1    3000    4000    peak_2    150    +
```

### Annotation Files

Place genome annotation files in the annotation directory or specify with `--custom`.
Supported genomes: HG18, HG19, HG38, MM8, MM9, MM10

## Output Files

ROSE generates several output files:

1. **`*_STITCHED.gff`**: Stitched enhancer regions
2. **`*_REGION_MAP.txt`**: Table with region coordinates and signal densities
3. **`*_SuperStitched.table.txt`**: Super-enhancer table
4. **`*_AllStitched.table.txt`**: All stitched enhancers
5. **`*_REGION_TO_GENE.txt`**: Enhancers mapped to genes
6. **`*_GENE_TO_REGION.txt`**: Genes mapped to enhancers
7. **PDF plots**: Visualization of super-enhancers

## Examples

### Basic Analysis

```bash
rose2 main \
    -i H3K27ac_peaks.gff \
    -r H3K27ac.bam \
    -g HG19 \
    -o rose_output/
```

### With Control and Custom Stitching

```bash
rose2 main \
    -i H3K27ac_peaks.gff \
    -r H3K27ac.bam \
    -c Input.bam \
    -g HG19 \
    -s 12500 \
    -t 2500 \
    -o rose_output/
```

### Using Custom Genome

```bash
rose2 main \
    -i peaks.gff \
    -r sample.bam \
    --custom /path/to/custom_genome.ucsc \
    -o rose_output/
```

## Python API

You can also use ROSE2 as a Python library:

```python
from rose2 import utils

# Create a locus
locus = utils.Locus("chr1", 1000, 2000, "+", "peak_1")

# Work with BAM files
bam = utils.Bam("sample.bam")
reads = bam.get_reads_locus(locus, sense="both")

# Parse GFF files
collection = utils.gff_to_locus_collection("peaks.gff")
stitched = collection.stitch_collection(stitch_window=12500)
```

## Credits

### Original Algorithm
- **Richard Young Lab**, Whitehead Institute for Biomedical Research
- **Contact**: youngcomputation@wi.mit.edu

### Python 3 Port
- **St. Jude Children's Research Hospital**, Abra Lab

### Modernization & PyPI Package
- **Ming (Tommy) Tang**
- **Contact**: tangming2005@gmail.com

## Citation

If you use ROSE2 in your research, please cite:

> Whyte et al. (2013) "Master Transcription Factors and Mediator Establish Super-Enhancers at Key Cell Identity Genes"
> *Cell* 153(2):307-319
> [https://doi.org/10.1016/j.cell.2013.03.035](https://doi.org/10.1016/j.cell.2013.03.035)

## License

Apache License 2.0 - see [LICENSE.txt](LICENSE.txt) for details

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions or issues:
- Open an issue on [GitHub](https://github.com/stjude/ROSE2/issues)
- Contact: tangming2005@gmail.com

## Changelog

### Version 2.0.0 (2025)
- Complete modernization to Python 3.8+ with type hints
- Replaced deprecated `optparse` with `argparse`
- Replaced `os.system()` with `subprocess` for better security
- Added proper logging throughout
- PyPI-installable package
- Improved error handling and validation
- Added comprehensive documentation
- Added test suite

### Previous Versions
- Python 3 port by St. Jude Children's Research Hospital
- Original Python 2 version by Young Lab

## Links

- **PyPI**: https://pypi.org/project/rose2/
- **GitHub**: https://github.com/stjude/ROSE2
- **Original Source**: https://bitbucket.org/young_computation/rose/
- **Paper**: https://pmc.ncbi.nlm.nih.gov/articles/PMC3841062/
