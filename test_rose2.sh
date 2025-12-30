#!/bin/bash
# ROSE2 Test Script
# This script shows how to run ROSE2 with real data

set -e  # Exit on error

echo "=== ROSE2 Test Script ==="
echo ""

# Check if ROSE2 is installed
if ! command -v rose2 &> /dev/null; then
    echo "ERROR: ROSE2 not found. Please install first:"
    echo "  pip install -e ."
    exit 1
fi

echo "✓ ROSE2 installed: $(rose2 --version)"
echo ""

# Check dependencies
echo "Checking dependencies..."

if ! command -v samtools &> /dev/null; then
    echo "WARNING: samtools not found. Install with:"
    echo "  brew install samtools  # macOS"
    echo "  sudo apt-get install samtools  # Ubuntu"
fi

if ! command -v Rscript &> /dev/null; then
    echo "WARNING: R not found. Install with:"
    echo "  brew install r  # macOS"
    echo "  sudo apt-get install r-base  # Ubuntu"
fi

if ! command -v bedtools &> /dev/null; then
    echo "WARNING: bedtools not found. Install with:"
    echo "  brew install bedtools  # macOS"
    echo "  sudo apt-get install bedtools  # Ubuntu"
fi

echo ""
echo "=== Example ROSE2 Commands ==="
echo ""

# Example 1: Show help
echo "1. Getting help:"
echo "   rose2 --help"
echo "   rose2 main --help"
echo ""

# Example 2: Basic run
echo "2. Basic run (requires your data):"
echo "   rose2 main -i peaks.gff -r sample.bam -g HG19 -o output/"
echo ""

# Example 3: With control
echo "3. With control BAM:"
echo "   rose2 main -i peaks.gff -r sample.bam -c control.bam -g HG19 -o output/"
echo ""

# Example 4: Custom parameters
echo "4. Custom parameters:"
echo "   rose2 main \\"
echo "       -i H3K27ac_peaks.gff \\"
echo "       -r H3K27ac.bam \\"
echo "       -c Input.bam \\"
echo "       -g HG19 \\"
echo "       -s 12500 \\"
echo "       -t 2500 \\"
echo "       -o rose_output/"
echo ""

# Example 5: Download test data
echo "5. Download example ChIP-seq data:"
echo ""
echo "   # From ENCODE (example):"
echo "   wget https://www.encodeproject.org/files/ENCFF000XYZ/@@download/ENCFF000XYZ.bam"
echo "   wget https://www.encodeproject.org/files/ENCFF000ABC/@@download/ENCFF000ABC.bed"
echo ""
echo "   # Convert BED to GFF:"
echo "   awk 'BEGIN{OFS=\"\\t\"} {print \$1,\$4,\".\",\$2,\$3,\$5,\$6,\".\",\$4}' peaks.bed > peaks.gff"
echo ""
echo "   # Sort and index BAM:"
echo "   samtools sort input.bam -o sorted.bam"
echo "   samtools index sorted.bam"
echo ""

echo "=== Ready to run ROSE2! ==="
echo ""
echo "Next steps:"
echo "1. Prepare your input files (peaks + BAM)"
echo "2. Make sure BAM is sorted and indexed"
echo "3. Run: rose2 main -i <peaks> -r <bam> -g <genome> -o <output>"
echo ""
