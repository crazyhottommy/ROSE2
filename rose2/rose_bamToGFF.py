#!/usr/bin/env python3
"""
BAM to GFF Mapper

Maps reads from a BAM file to genomic regions defined in a GFF file and
calculates read density across those regions.

Original: Young Lab, Whitehead Institute
Python 3 Port: St. Jude Children's Research Hospital
Modernization: Ming Tang
"""

import argparse
import logging
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, DefaultDict, List, Optional, Union

from rose2 import utils

logger = logging.getLogger(__name__)


def map_bam_to_gff(
    bam_file: str,
    gff: Union[str, List[List[str]]],
    sense: str = "both",
    extension: int = 200,
    floor: int = 0,
    rpm: bool = False,
    matrix: Optional[int] = None,
) -> List[List[Any]]:
    """Map reads from a BAM file to GFF regions.

    Args:
        bam_file: Path to sorted and indexed BAM file
        gff: Either path to GFF file or parsed GFF table
        sense: Strand to consider ('+', '-', or 'both')
        extension: Extend reads by this many bp (default: 200)
        floor: Minimum read threshold to count towards density
        rpm: If True, normalize to reads per million
        matrix: If specified, output fixed-bin-number matrix

    Returns:
        New GFF table with density information
    """
    floor = int(floor)

    # Initialize BAM object
    bam = utils.Bam(bam_file)

    # Output table
    new_gff: List[List[Any]] = []

    # Calculate normalization factor
    if rpm:
        mmr = round(float(bam.get_total_reads("mapped")) / 1000000, 4)
    else:
        mmr = 1

    logger.info(f"Using MMR value of {mmr}")

    # Check if BAM uses 'chr' prefix
    has_chr_flag = utils.check_chr_status(bam_file)
    if has_chr_flag:
        logger.info("BAM file has 'chr' prefix")
    else:
        logger.info("BAM file does not have 'chr' prefix")

    # Parse GFF if it's a file path
    if isinstance(gff, str):
        gff = utils.parse_table(gff, "\t")

    # Set up matrix table header
    if matrix:
        new_gff.append(
            ["GENE_ID", "locusLine"]
            + [f"bin_{n}_{bam_file.split('/')[-1]}" for n in range(1, int(matrix) + 1)]
        )

    # Process each GFF line
    ticker = 0
    logger.info("Processing GFF regions")
    for line in gff:
        line = line[0:9]
        if ticker % 100 == 0:
            logger.info(f"Processed {ticker} regions")
        ticker += 1

        # Remove 'chr' prefix if BAM doesn't have it
        if not has_chr_flag:
            line[0] = re.sub(r"chr", r"", line[0])

        gff_locus = utils.Locus(line[0], int(line[3]), int(line[4]), line[6], line[1])
        search_locus = utils.make_search_locus(gff_locus, int(extension), int(extension))

        # Get reads for this locus
        reads = bam.get_reads_locus(search_locus, "both", False, "none")

        # Extend reads
        extended_reads = []
        for locus in reads:
            if locus.sense() == "+" or locus.sense() == ".":
                locus = utils.Locus(
                    locus.chr(),
                    locus.start(),
                    locus.end() + extension,
                    locus.sense(),
                    locus.ID(),
                )
            if locus.sense() == "-":
                locus = utils.Locus(
                    locus.chr(),
                    locus.start() - extension,
                    locus.end(),
                    locus.sense(),
                    locus.ID(),
                )
            extended_reads.append(locus)

        # Separate sense and antisense reads
        if gff_locus.sense() == "+" or gff_locus.sense() == ".":
            sense_reads = [x for x in extended_reads if x.sense() == "+" or x.sense() == "."]
            anti_reads = [x for x in extended_reads if x.sense() == "-"]
        else:
            sense_reads = [x for x in extended_reads if x.sense() == "-" or x.sense() == "."]
            anti_reads = [x for x in extended_reads if x.sense() == "+"]

        # Create read density hashes
        sense_hash: DefaultDict[int, int] = defaultdict(int)
        anti_hash: DefaultDict[int, int] = defaultdict(int)

        # Fill in read hashes
        if sense in ["+", "both", "."]:
            for read in sense_reads:
                for x in range(read.start(), read.end() + 1):
                    sense_hash[x] += 1

        if sense in ["-", "both", "."]:
            for read in anti_reads:
                for x in range(read.start(), read.end() + 1):
                    anti_hash[x] += 1

        # Apply flooring and coordinate filtering
        keys = utils.uniquify(list(sense_hash.keys()) + list(anti_hash.keys()))
        if floor > 0:
            keys = [x for x in keys if (sense_hash[x] + anti_hash[x]) > floor]

        # Filter for coordinates within GFF locus
        keys = [x for x in keys if gff_locus.start() < x < gff_locus.end()]

        # Set up output line
        if not has_chr_flag:
            cluster_line = [gff_locus.ID(), "chr" + str(gff_locus)]
        else:
            cluster_line = [gff_locus.ID(), str(gff_locus)]

        # Calculate bin density
        if matrix:
            bin_size = (gff_locus.len() - 1) / int(matrix)
            n_bins = int(matrix)

            if bin_size == 0:
                cluster_line += ["NA"] * int(matrix)
                new_gff.append(cluster_line)
                continue

            n = 0
            if gff_locus.sense() in ["+", ".", "both"]:
                i = gff_locus.start()
                while n < n_bins:
                    n += 1
                    bin_keys = [x for x in keys if i < x < i + bin_size]
                    bin_den = float(sum([sense_hash[x] + anti_hash[x] for x in bin_keys])) / bin_size
                    cluster_line.append(round(bin_den / mmr, 4))
                    i = i + bin_size
            else:
                i = gff_locus.end()
                while n < n_bins:
                    n += 1
                    bin_keys = [x for x in keys if i - bin_size < x < i]
                    bin_den = float(sum([sense_hash[x] + anti_hash[x] for x in bin_keys])) / bin_size
                    cluster_line.append(round(bin_den / mmr, 4))
                    i = i - bin_size

        new_gff.append(cluster_line)

    return new_gff


def main(argv: Optional[List[str]] = None) -> int:
    """Main execution function.

    Args:
        argv: Command line arguments (for testing)

    Returns:
        Exit code
    """
    # Set up logging
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    parser = argparse.ArgumentParser(
        description="Map BAM reads to GFF regions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument("-b", "--bam", required=True, help="Sorted and indexed BAM file")
    required.add_argument(
        "-i", "--input", required=True, help="Input .gff or enriched region file"
    )

    # Optional arguments
    parser.add_argument("-o", "--output", help="Output filename (default: input + .mapped)")
    parser.add_argument(
        "-s",
        "--sense",
        choices=["+", "-", ".", "both"],
        default="both",
        help="Map to '+', '-', or 'both' strands (default: both)",
    )
    parser.add_argument(
        "-f",
        "--floor",
        type=int,
        default=0,
        help="Read floor threshold necessary to count towards density (default: 0)",
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=int,
        default=200,
        help="Extend reads by N bp (default: 200)",
    )
    parser.add_argument(
        "-r",
        "--rpm",
        action="store_true",
        default=False,
        help="Normalize density to reads per million",
    )
    parser.add_argument(
        "-m",
        "--matrix",
        type=int,
        help="Output variable bin sized matrix with specified number of bins",
    )

    args = parser.parse_args(argv)

    # Validate BAM file
    bam_path = Path(args.bam).resolve()
    if not bam_path.exists():
        logger.error(f"BAM file not found: {args.bam}")
        return 1

    # Check for BAI index
    bai_path = Path(str(bam_path) + ".bai")
    bam_dir = bam_path.parent
    bam_name = bam_path.stem

    has_bai = False
    for file_path in bam_dir.iterdir():
        if file_path.name.startswith(bam_name) and file_path.suffix == ".bai":
            has_bai = True
            break

    if not has_bai:
        logger.error("No associated .bai index file found. BAM must be sorted and indexed.")
        return 1

    # Validate input file
    if not Path(args.input).exists():
        logger.error(f"Input file not found: {args.input}")
        return 1

    # Set output path
    if args.output:
        output = args.output
    else:
        output = os.getcwd() + "/" + Path(args.input).name + ".mapped"

    # Run mapping
    logger.info("Mapping BAM to GFF and creating matrix with fixed bin number")
    new_gff = map_bam_to_gff(
        args.bam,
        args.input,
        args.sense,
        args.extension,
        args.floor,
        args.rpm,
        args.matrix,
    )

    # Write output
    utils.unparse_table(new_gff, output, "\t")
    logger.info(f"Output written to {output}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
