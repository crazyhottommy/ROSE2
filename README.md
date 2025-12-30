# ROSE2: Rank Ordering of Super-Enhancers

[![PyPI version](https://badge.fury.io/py/rose2.svg)](https://badge.fury.io/py/rose2)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

A fast, modern tool for identifying super-enhancers and their target genes from ChIP-seq data.

---

## What Problem Does ROSE2 Solve?

### The Biological Challenge

Cells maintain their identity and function through **gene regulatory networks** controlled by enhancers—DNA elements that activate gene expression. Among these, **super-enhancers** are especially important:

- **Cell identity genes**: Super-enhancers drive the expression of genes that define cell type (e.g., MYC in cancer cells, OCT4 in stem cells)
- **Disease-associated genes**: Many disease genes are controlled by super-enhancers, making them therapeutic targets
- **Master regulators**: Super-enhancers regulate transcription factors that control entire gene programs

However, identifying which genomic regions are super-enhancers is challenging because:
1. Enhancers are scattered across the genome, often far from their target genes
2. Multiple enhancer elements often work together as clusters
3. Standard analysis treats each enhancer separately, missing the biological reality of collaborative regulation

### Why Super-Enhancers Matter

Super-enhancers are:

🎯 **Therapeutic targets**: BET inhibitors and other drugs target super-enhancer components. 
🔬 **Biomarkers**: Super-enhancer landscapes distinguish cell types and disease states. 
🧬 **Mechanistic insights**: Reveal how transcription factors and cofactors coordinate gene expression. 
💊 **Drug development**: Understanding super-enhancers helps design selective therapies. 

**Example applications:**
- **Cancer biology**: Identify oncogene-driving super-enhancers (MYC, NOTCH1, RUNX1)
- **Stem cell research**: Map pluripotency-associated super-enhancers (OCT4, NANOG, SOX2)
- **Immunology**: Discover immune cell identity super-enhancers (CD4, IL2RA, FOXP3)
- **Development**: Track super-enhancer dynamics during differentiation

---

## How ROSE2 Solves It

ROSE2 identifies super-enhancers through a **biologically-motivated algorithm**:

### 1. **Stitch nearby enhancers** (the biological reality)
   - Combines enhancer peaks within 12.5kb (default)
   - Reflects how enhancers physically cluster in 3D chromatin space
   - Creates "stitched regions" representing functional enhancer units

### 2. **Rank by regulatory activity** (quantitative measurement)
   - Calculates ChIP-seq signal density for each stitched region
   - Accounts for input/control background
   - Normalizes for region size

### 3. **Identify super-enhancers** (data-driven threshold)
   - Plots signal distribution (hockey stick curve)
   - Super-enhancers are the inflection point: exceptionally high signal
   - Typically the top 5-10% of stitched enhancers

### 4. **Map to target genes** (biological function)
   - Links enhancers to nearby genes (overlapping, proximal, or closest)
   - Enables functional interpretation
   - Identifies putative regulatory targets

**Input**: ChIP-seq data (BAM files) + peak calls (BED/GFF) for enhancer-associated factors
- **Commonly used**: H3K27ac (active enhancer mark), BRD4, MED1, p300
- Works with any ChIP-seq data for enhancer-binding factors

**Output**: Ranked enhancer list, super-enhancer calls, gene assignments, plots

📖 **For algorithmic details, see [TECHNICAL_NOTES.md](TECHNICAL_NOTES.md)**

---

## Installation

### Quick Install (Recommended)

```bash
pip install pyrose2
```

**Note:** The package is installed as `pyrose2`, but the command-line tool is still `rose2`:
- Install: `pip install pyrose2`
- Use: `rose2 main -i peaks.bed ...`

### System Requirements

**Required:**
- Python 3.8 or higher
- [samtools](http://www.htslib.org/) - for BAM file processing
- [R](https://www.r-project.org/) (≥ 3.4) - for plotting and statistics

**Optional but recommended:**
- [bedtools](https://bedtools.readthedocs.io/) (≥ 2) - for format conversions

**Install dependencies:**

Ubuntu/Debian:
```bash
sudo apt-get install samtools bedtools r-base
```

macOS (Homebrew):
```bash
brew install samtools bedtools r
```

Conda:
```bash
conda install -c bioconda samtools bedtools r-base
```

---

## Quick Start

### Example: Identify H3K27ac Super-Enhancers

```bash
# Basic usage
rose2 main \
  -i H3K27ac_peaks.bed \
  -r H3K27ac.bam \
  -g HG38 \
  -o results/

# With input control
rose2 main \
  -i H3K27ac_peaks.bed \
  -r H3K27ac.bam \
  -c Input.bam \
  -g HG38 \
  -o results/

# Custom stitching distance (e.g., 20kb)
rose2 main \
  -i BRD4_peaks.bed \
  -r BRD4.bam \
  -g MM10 \
  -s 20000 \
  -o results/
```

### Input Files

**1. Peak calls** (BED, narrowPeak, or GFF format):
```
chr1    1000    2000    peak1    100    .
chr1    5000    6000    peak2    150    .
```

**2. Aligned reads** (sorted, indexed BAM file):
```bash
samtools sort input.bam -o sorted.bam
samtools index sorted.bam
```

**3. Genome annotation**: Built-in support for HG18, HG19, HG38, MM8, MM9, MM10

---

## Output Files

ROSE2 generates several files in your output directory:

### Core Results

**1. Super-enhancer table** (`*_SuperStitched.table.txt`)
- List of identified super-enhancers
- Genomic coordinates, signal strength, nearby genes
- **Key for downstream analysis**

**2. All stitched enhancers** (`*_AllStitched.table.txt`)
- Complete ranked list of all stitched regions
- Allows custom thresholding

**3. Hockey stick plot** (`*_Plot_points.png`)
- Visual identification of super-enhancers
- Shows signal distribution and threshold

### Gene Mapping Results

**4. Enhancer-to-gene mapping** (`*_REGION_TO_GENE.txt`)
- Each row = one enhancer
- Shows overlapping, proximal, and closest genes

**5. Gene-to-enhancer mapping** (`*_GENE_TO_REGION.txt`)
- Each row = one gene
- Lists all associated enhancers

**6. Signal with genes** (`*.table_withGENES.txt`)
- Combines enhancer signal with gene assignments
- **Ready for gene set enrichment analysis**

---

## Typical Workflows

### 1. Basic Super-Enhancer Discovery

**Goal**: Identify super-enhancers in your cell type

```bash
# Run ROSE2
rose2 main -i peaks.bed -r H3K27ac.bam -g HG38 -o results/

# Analyze results
# - Check hockey stick plot for super-enhancer threshold
# - Review super-enhancer gene list for cell identity genes
# - Compare to other cell types or conditions
```

### 2. Compare Conditions (e.g., Cancer vs Normal)

```bash
# Run on both conditions
rose2 main -i cancer_peaks.bed -r cancer.bam -g HG38 -o cancer_SE/
rose2 main -i normal_peaks.bed -r normal.bam -g HG38 -o normal_SE/

# Compare outputs
# - Identify cancer-specific super-enhancers
# - Check for oncogene associations (MYC, NOTCH1, etc.)
# - Analyze gained/lost super-enhancers
```

### 3. Time Course (e.g., Differentiation)

```bash
# Run on each timepoint
for day in day0 day2 day4 day7; do
  rose2 main -i ${day}_peaks.bed -r ${day}.bam -g MM10 -o ${day}_SE/
done

# Track super-enhancer dynamics
# - Identify stage-specific super-enhancers
# - Map to developmental regulators
# - Build regulatory trajectories
```

### 4. Multi-Factor Analysis (e.g., BRD4 + H3K27ac)

```bash
# Run with different ChIP targets
rose2 main -i peaks.bed -r H3K27ac.bam -g HG38 -o H3K27ac_SE/
rose2 main -i peaks.bed -r BRD4.bam -g HG38 -o BRD4_SE/
rose2 main -i peaks.bed -r MED1.bam -g HG38 -o MED1_SE/

# Compare factor occupancy
# - Identify co-occupied super-enhancers
# - Assess factor dependency
# - Predict drug sensitivity
```

---

## Interpreting Results

### Understanding the Hockey Stick Plot

The plot shows cumulative ChIP-seq signal vs. rank:
- **X-axis**: Enhancers ranked by signal (1 = strongest)
- **Y-axis**: Cumulative ChIP-seq signal
- **Inflection point**: Where the curve "bends" sharply upward
- **Super-enhancers**: Regions above the inflection point (typically top 5-10%)

![Example Hockey Stick Plot](brd4_super.png)
*Example: BRD4 super-enhancer analysis showing the characteristic "hockey stick" curve. The red line indicates the inflection point separating typical enhancers (left) from super-enhancers (right).*

**What to look for:**
- Clear separation between typical and super-enhancers
- Super-enhancers should be 5-10× stronger than typical enhancers
- Too many super-enhancers? Increase stitching distance
- Too few? Decrease stitching distance or check data quality

### Gene Assignment Strategy

ROSE2 assigns genes using a proximity-based hierarchy:

1. **Overlapping genes**: Gene body overlaps the enhancer (distance = 0)
2. **Proximal genes**: TSS within 50kb of enhancer boundary
3. **Closest gene**: Nearest TSS if no overlapping/proximal genes found

**Biological interpretation:**
- Overlapping = intragenic enhancer (common for large genes)
- Proximal = typical enhancer-promoter distance
- Distal = long-range regulation (validate with Hi-C/3C data)

### Validating Super-Enhancers

**Computational validation:**
- Check for cell identity genes (expected master regulators)
- Compare to published super-enhancers in similar cell types
- Verify enrichment at known regulatory loci

**Experimental validation:**
- CRISPR deletion of super-enhancer regions
- BET inhibitor treatment (should reduce super-enhancer activity)
- 3C/Hi-C to confirm enhancer-promoter interactions
- RNA-seq after enhancer perturbation

---

## Advanced Usage

### Custom Parameters

```bash
rose2 main \
  -i peaks.bed \
  -r sample.bam \
  -c input.bam \
  -g HG38 \
  -s 12500 \          # Stitching distance (default: 12.5kb)
  -t 2500 \           # TSS exclusion zone (default: 2.5kb)
  -o output/ \
  --mask blacklist.bed # Exclude problematic regions
```

**Parameter guidelines:**
- **Stitching distance (-s)**: Larger = fewer, bigger super-enhancers
  - 12.5kb (default): Standard for most analyses
  - 20kb: For very dense enhancer regions
  - 5kb: For sparse enhancer landscapes

- **TSS exclusion (-t)**: Exclude promoter-proximal regions
  - 2.5kb (default): Removes promoters while keeping enhancers
  - 0: Include all regions (not recommended)
  - 5kb: More stringent enhancer-only analysis

### Custom Genome Annotation

```bash
# Use your own gene annotation
rose2 main \
  -i peaks.bed \
  -r sample.bam \
  --custom my_annotation.ucsc \
  -o output/
```

**Annotation format** (UCSC refGene format):
```
585    NR_046018    chr1    +    11873    14409    14409    14409    3    11873,12612,13220,    12227,12721,14409,    0    DDX11L1    unk    unk    -1,-1,-1,
```

### Gene List Filtering

```bash
# Map only to specific genes of interest
rose2-geneMapper \
  -i SuperStitched.table.txt \
  -g HG38 \
  -l my_genes.txt \     # One gene per line
  -o filtered_mapping/
```

**Use cases:**
- Focus on known oncogenes/tumor suppressors
- Analyze specific pathways (e.g., immune genes)
- Validate predictions for candidate genes

---

## Performance

ROSE2 v2.0 is **dramatically faster** than previous versions:

| Dataset | Previous Version | ROSE2 v2.0 | Speedup |
|---------|-----------------|------------|---------|
| 22K enhancers, HG38 | ~9 hours | ~12 minutes | **45×** |
| Gene mapping | 6 hours | 12 seconds | **1,772×** |
| Memory usage | 5 GB | 500 MB | **90% reduction** |

**Real-world example (measured on MOLM-14 BRD4 ChIP-seq):**
- **Runtime**: 1.7 minutes
- **Peak memory**: 620 MB
- **Input**: 22,138 narrowPeak regions (MACS3 output)
- **Genome**: HG38 with 58,037 RefSeq transcripts

**Optimizations:**
- Interval-based coverage calculation (500× faster)
- Smart gene search algorithm (eliminates billions of redundant checks)
- Streaming data processing (90% memory reduction)
- See [CHANGELOG.md](CHANGELOG.md) for details

---

## Troubleshooting

### Common Issues

**1. "samtools not found"**
```bash
# Install samtools
conda install -c bioconda samtools
# OR
brew install samtools  # macOS
sudo apt-get install samtools  # Ubuntu
```

**2. "BAM file not indexed"**
```bash
samtools index your_file.bam
# Creates your_file.bam.bai
```

**3. "No super-enhancers found"**
- Check data quality (coverage, signal-to-noise)
- Verify peak calls are reasonable
- Try adjusting stitching distance
- Ensure you're using enhancer marks (H3K27ac, not H3K4me3)

**4. "Wrong chromosome naming (chr1 vs 1)"**
- ROSE2 automatically detects and handles this
- If issues persist, check BAM header: `samtools view -H your.bam`

**5. "Out of memory"**
- ROSE2 v2.0 uses 90% less memory than before
- If still issues, process smaller chromosome regions separately

---

## Citation

If you use ROSE2 in your research, please cite:

**Original ROSE algorithm:**
> Whyte WA, Orlando DA, Hnisz D, et al. Master transcription factors and mediator establish super-enhancers at key cell identity genes. *Cell*. 2013;153(2):307-319. doi:10.1016/j.cell.2013.03.035

**ROSE2 modernization and optimization:**
> Tang M (2025). ROSE2: High-performance super-enhancer identification for Python 3. https://github.com/crazyhottommy/ROSE2

---

## Credits

**Original Algorithm**: Richard Young Lab, Whitehead Institute.   
**Python 3 Port**: St. Jude Children's Research Hospital, Abra Lab.   
**Modernization & Optimization**: Ming (Tommy) Tang.   

- Modern Python packaging and PyPI distribution
- 1,700× performance improvements
- 90% memory reduction
- Comprehensive documentation and testing

---

## Related Tools

- **HOMER** - Motif analysis and peak annotation
- **GREAT** - Functional enrichment of regulatory regions
- **deepTools** - ChIP-seq quality control and visualization
- **ChromHMM** - Chromatin state discovery and characterization

---

## Getting Help

- **Documentation**: [Full guide](TECHNICAL_NOTES.md)
- **Issues**: [GitHub Issues](https://github.com/crazyhottommy/ROSE2/issues)
- **Email**: tangming2005@gmail.com

---

## License

Apache License 2.0 - See [LICENSE.txt](LICENSE.txt) for details.

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history and detailed changes.
