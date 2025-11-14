
# ROSE: Rank Ordering of Super-Enhancers

ROSE (Rank Ordering of Super-Enhancers) is a computational pipeline for the identification and analysis of super-enhancers from ChIP-seq data. This repository provides a Python 3-compatible, location-independent implementation, with additional support for Docker-based deployment.

## Overview

ROSE enables researchers to:
- Identify super-enhancers from ChIP-seq data
- Rank enhancer regions by signal density
- Map stitched enhancers to genes

## Key Changes

- Compatible with Python 3
- Executable independent of installation directory
- Available as a Docker image: `ghcr.io/stjude/abralab/rose:latest`

## Installation & Setup

Clone the repository and set up the environment variables:

```bash
ROSEPATH=/path/to/ROSE
export PYTHONPATH="$PYTHONPATH:$ROSEPATH/lib"
export PATH="$PATH:$ROSEPATH/bin"
```

## Usage

Run the main program with the following command:

```bash
ROSE_main.py -g [GENOME] -i [INPUT_REGION_GFF] -r [RANKBY_BAM_FILE] -o [OUTPUT_FOLDER] [OPTIONAL_FLAGS]
```

## Dependencies

- Python 3
- [samtools](http://www.htslib.org/)
- R (version > 3.4)
- bedtools (version > 2)

## Input Requirements

All input files must reside in a single directory.

### Annotation File (`-g [GENOME]`)
- Must be in UCSC table track format ([UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables))
- Filename: `[GENOME]_refseq.ucsc` (e.g., `hg19_refseq.ucsc`)
- Place in the `annotation/` folder within the input directory

### BAM File (`-r [RANKBY_BAM_FILE]`)
- Chromosome IDs must start with `chr`
- Files must be coordinate sorted and indexed using SAMtools

### Peak File (Constituent Enhancers, `-i [INPUT_REGION_GFF]`)
- Format: GFF
- Required columns:
	- Column 1: Chromosome (chr#)
	- Column 2: Unique ID for each constituent enhancer region
	- Column 4: Start position
	- Column 5: End position
	- Column 7: Strand (`+`, `-`, or `.`)
	- Column 9: Unique ID for each constituent enhancer region
- **Note:** If values in columns 2 and 9 differ, column 2 will be used.

## Directory Structure

```
├── LICENSE.txt
├── README.md
├── bin/
│   ├── ROSE_bamToGFF.py      # Calculates density of .bam reads in .gff regions
│   ├── ROSE_callSuper.R      # Ranks regions by their densities, creates cutoff
│   ├── ROSE_geneMapper.py    # Assigns stitched enhancers to genes
│   └── ROSE_main.py          # Main program
└── lib/
	└── ROSE_utils.py         # Utility methods
```

## References

- Original source: [Bitbucket - young_computation/rose](https://bitbucket.org/young_computation/rose/src/master/)

For questions or support, please contact the repository maintainers.
