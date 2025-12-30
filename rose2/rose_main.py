#!/usr/bin/env python3
"""
ROSE Main Pipeline

Stitches together regions to form enhancers, maps read density to stitched regions,
and ranks enhancers by read density to discover super-enhancers.

Original: Young Lab, Whitehead Institute
Python 3 Port: St. Jude Children's Research Hospital
Modernization: Ming Tang
"""

import argparse
import logging
import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, List, Optional, Tuple

from rose2 import utils

logger = logging.getLogger(__name__)


def region_stitching(
    input_gff: str,
    stitch_window: int,
    tss_window: int,
    annot_file: str,
    remove_tss: bool = True,
) -> Tuple[utils.LocusCollection, List[List[Any]]]:
    """Stitch together genomic regions to form enhancers.

    Args:
        input_gff: Path to input GFF file
        stitch_window: Maximum distance for stitching regions
        tss_window: Distance from TSS to exclude
        annot_file: Path to annotation file
        remove_tss: If True, remove TSS-overlapping regions

    Returns:
        Tuple of (stitched LocusCollection, debug output)
    """
    logger.info("PERFORMING REGION STITCHING")

    bound_collection = utils.gff_to_locus_collection(input_gff)
    debug_output: List[List[Any]] = []

    if remove_tss:
        # Create TSS loci for filtering
        start_dict = utils.make_start_dict(annot_file)
        tss_loci = []
        for gene_id in list(start_dict.keys()):
            tss_loci.append(utils.make_tss_locus(gene_id, start_dict, tss_window, tss_window))

        tss_collection = utils.LocusCollection(tss_loci, 50)
        bound_loci = bound_collection.get_loci()

        # Remove regions contained by TSS exclusion zone
        remove_ticker = 0
        for locus in bound_loci:
            if len(tss_collection.get_containers(locus, "both")) > 0:
                bound_collection.remove(locus)
                debug_output.append([str(locus), locus.ID(), "CONTAINED"])
                remove_ticker += 1

        logger.info(f"REMOVED {remove_ticker} LOCI BECAUSE THEY WERE CONTAINED BY A TSS")

    # Stitch the collection
    stitched_collection = bound_collection.stitch_collection(stitch_window, "both")

    if remove_tss:
        # Replace stitched regions overlapping multiple genes with original loci
        fixed_loci = []
        tss_loci = []
        start_dict = utils.make_start_dict(annot_file)
        for gene_id in list(start_dict.keys()):
            tss_loci.append(utils.make_tss_locus(gene_id, start_dict, 50, 50))

        tss_collection = utils.LocusCollection(tss_loci, 50)
        remove_ticker = 0
        original_ticker = 0

        for stitched_locus in stitched_collection.get_loci():
            overlapping_tss_loci = tss_collection.get_overlap(stitched_locus, "both")
            tss_names = [
                start_dict[tss_locus.ID()]["name"] for tss_locus in overlapping_tss_loci
            ]
            tss_names = utils.uniquify(tss_names)

            if len(tss_names) > 2:
                original_loci = bound_collection.get_overlap(stitched_locus, "both")
                original_ticker += len(original_loci)
                fixed_loci += original_loci
                debug_output.append([str(stitched_locus), stitched_locus.ID(), "MULTIPLE_TSS"])
                remove_ticker += 1
            else:
                fixed_loci.append(stitched_locus)

        logger.info(f"REMOVED {remove_ticker} STITCHED LOCI BECAUSE THEY OVERLAPPED MULTIPLE TSSs")
        logger.info(f"ADDED BACK {original_ticker} ORIGINAL LOCI")
        fixed_collection = utils.LocusCollection(fixed_loci, 50)
        return fixed_collection, debug_output
    else:
        return stitched_collection, debug_output


def map_collection(
    stitched_collection: utils.LocusCollection,
    reference_collection: utils.LocusCollection,
    bam_file_list: List[str],
    mapped_folder: str,
    output: str,
    ref_name: str,
) -> None:
    """Map factor density to stitched regions and create output table.

    Args:
        stitched_collection: Collection of stitched loci
        reference_collection: Collection of reference regions
        bam_file_list: List of BAM file paths
        mapped_folder: Folder containing mapped GFF files
        output: Output file path
        ref_name: Reference name for file naming
    """
    logger.info("FORMATTING TABLE")
    loci = stitched_collection.get_loci()

    # Remove chrY loci
    loci = [locus for locus in loci if locus.chr() != "chrY"]

    locus_table: List[List[Any]] = [
        ["REGION_ID", "CHROM", "START", "STOP", "NUM_LOCI", "CONSTITUENT_SIZE"]
    ]

    loci_len_list = [locus.len() for locus in loci]
    len_order = utils.order(loci_len_list, decreasing=True)

    ticker = 0
    for i in len_order:
        ticker += 1
        if ticker % 1000 == 0:
            logger.info(f"Processed {ticker} loci")

        locus = loci[i]

        # Get size of enriched regions within stitched locus
        ref_enrich_size = 0
        ref_overlapping_loci = reference_collection.get_overlap(locus, "both")
        for ref_locus in ref_overlapping_loci:
            ref_enrich_size += ref_locus.len()

        try:
            stitch_count = int(locus.ID().split("_")[0])
        except ValueError:
            stitch_count = 1

        locus_table.append(
            [locus.ID(), locus.chr(), locus.start(), locus.end(), stitch_count, ref_enrich_size]
        )

    logger.info("GETTING MAPPED DATA")
    for bam_file in bam_file_list:
        bam_file_name = bam_file.split("/")[-1]
        logger.info(f"GETTING MAPPING DATA FOR {bam_file}")

        mapped_gff_path = f"{mapped_folder}{ref_name}_{bam_file_name}_MAPPED.gff"
        logger.info(f"OPENING {mapped_gff_path}")

        mapped_gff = utils.parse_table(mapped_gff_path, "\t")

        signal_dict: DefaultDict[str, float] = defaultdict(float)
        logger.info(f"MAKING SIGNAL DICT FOR {bam_file}")

        mapped_loci = []
        for line in mapped_gff[1:]:
            chrom = line[1].split("(")[0]
            start = int(line[1].split(":")[-1].split("-")[0])
            end = int(line[1].split(":")[-1].split("-")[1])
            mapped_loci.append(utils.Locus(chrom, start, end, ".", line[0]))

            try:
                signal_dict[line[0]] = float(line[2]) * (abs(end - start))
            except ValueError:
                logger.warning(f"WARNING NO SIGNAL FOR LINE: {line}")
                continue

        mapped_collection = utils.LocusCollection(mapped_loci, 500)
        locus_table[0].append(bam_file_name)

        for i in range(1, len(locus_table)):
            signal = 0.0
            line = locus_table[i]
            line_locus = utils.Locus(line[1], line[2], line[3], ".")
            overlapping_regions = mapped_collection.get_overlap(line_locus, sense="both")
            for region in overlapping_regions:
                signal += signal_dict[region.ID()]
            locus_table[i].append(signal)

    utils.unparse_table(locus_table, output, "\t")


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
        description="ROSE: Rank Ordering of Super-Enhancers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input file of binding sites (.gff, .bed, .narrowPeak, or .broadPeak)",
    )
    required.add_argument(
        "-r", "--rankby", required=True, help="BAM file to rank enhancers by"
    )
    required.add_argument("-o", "--out", required=True, help="Output folder")

    # Genome arguments (one required)
    genome_group = parser.add_mutually_exclusive_group(required=True)
    genome_group.add_argument(
        "-g",
        "--genome",
        choices=["MM9", "MM8", "MM10", "HG18", "HG19", "HG38"],
        help="Genome build",
    )
    genome_group.add_argument(
        "--custom", dest="custom_genome", help="Custom genome annotation refseq.ucsc file"
    )

    # Optional arguments
    parser.add_argument(
        "-b",
        "--bams",
        help="Comma-separated list of additional BAM files to map to",
    )
    parser.add_argument("-c", "--control", help="Control BAM file")
    parser.add_argument(
        "-s",
        "--stitch",
        type=int,
        default=12500,
        help="Maximum linking distance for stitching (default: 12500)",
    )
    parser.add_argument(
        "-t",
        "--tss",
        type=int,
        default=0,
        help="Distance from TSS to exclude. 0 = no TSS exclusion (default: 0)",
    )

    args = parser.parse_args(argv)

    # Create output folder
    out_folder = utils.format_folder(args.out, True)
    if not out_folder:
        logger.error(f"Could not create output folder: {args.out}")
        return 1

    # Set up folder structure
    gff_folder = utils.format_folder(os.path.join(out_folder, "gff"), True)
    mapped_folder = utils.format_folder(os.path.join(out_folder, "mappedGFF"), True)

    # Process input file
    if args.input.endswith(".bed") or args.input.endswith(".narrowPeak") or args.input.endswith(".broadPeak"):
        # Convert BED/narrowPeak/broadPeak to GFF
        input_gff_name = Path(args.input).stem
        input_gff_file = os.path.join(gff_folder, f"{input_gff_name}.gff")
        logger.info(f"Converting BED/Peak file to GFF format")
        utils.bed_to_gff(args.input, input_gff_file)
    elif args.input.endswith(".gff"):
        # Copy GFF to folder
        input_gff_file = args.input
        subprocess.run(["cp", input_gff_file, gff_folder], check=True)
    else:
        logger.warning("INPUT FILE DOES NOT END IN .gff, .bed, .narrowPeak, or .broadPeak. ASSUMING .gff FILE FORMAT")
        input_gff_file = args.input
        subprocess.run(["cp", input_gff_file, gff_folder], check=True)

    # Get list of BAM files to process (convert to absolute paths)
    bam_file_list = [str(Path(args.rankby).resolve())]
    if args.control:
        bam_file_list.append(str(Path(args.control).resolve()))
    if args.bams:
        bam_file_list += [str(Path(b).resolve()) for b in args.bams.split(",")]
        bam_file_list = utils.uniquify(bam_file_list)

    # Set stitching parameters
    stitch_window = args.stitch
    tss_window = args.tss
    remove_tss = tss_window != 0

    # Get input name
    logger.info(f"USING {input_gff_file} AS THE INPUT GFF")
    input_name = Path(input_gff_file).stem

    # Get annotation file (use package directory, not current directory)
    package_dir = Path(__file__).parent
    genome_dict = {
        "HG18": package_dir / "annotation" / "hg18_refseq.ucsc",
        "MM9": package_dir / "annotation" / "mm9_refseq.ucsc",
        "HG19": package_dir / "annotation" / "hg19_refseq.ucsc",
        "HG38": package_dir / "annotation" / "hg38_refseq.ucsc",
        "MM8": package_dir / "annotation" / "mm8_refseq.ucsc",
        "MM10": package_dir / "annotation" / "mm10_refseq.ucsc",
    }

    if args.custom_genome:
        annot_file = args.custom_genome
        logger.info(f"USING CUSTOM GENOME {args.custom_genome} AS THE GENOME FILE")
    else:
        genome = args.genome.upper()
        annot_file = str(genome_dict[genome])
        logger.info(f"USING {genome} AS THE GENOME")

    # Make start dict
    logger.info("MAKING START DICT")
    start_dict = utils.make_start_dict(annot_file)

    # Load reference collection
    logger.info("LOADING IN GFF REGIONS")
    reference_collection = utils.gff_to_locus_collection(input_gff_file)

    # Stitch regions
    logger.info("STITCHING REGIONS TOGETHER")
    stitched_collection, debug_output = region_stitching(
        input_gff_file, stitch_window, tss_window, annot_file, remove_tss
    )

    # Create stitched GFF
    logger.info("MAKING GFF FROM STITCHED COLLECTION")
    stitched_gff = utils.locus_collection_to_gff(stitched_collection)

    if not remove_tss:
        stitched_gff_file = os.path.join(
            gff_folder, f"{input_name}_{stitch_window // 1000}KB_STITCHED.gff"
        )
        stitched_gff_name = f"{input_name}_{stitch_window // 1000}KB_STITCHED"
    else:
        stitched_gff_file = os.path.join(
            gff_folder, f"{input_name}_{stitch_window // 1000}KB_STITCHED_TSS_DISTAL.gff"
        )
        stitched_gff_name = f"{input_name}_{stitch_window // 1000}KB_STITCHED_TSS_DISTAL"

    # Write stitched GFF
    logger.info(f"WRITING STITCHED GFF TO DISK AS {stitched_gff_file}")
    utils.unparse_table(stitched_gff, stitched_gff_file, "\t")

    # Set up output file
    output_file = os.path.join(out_folder, f"{stitched_gff_name}_REGION_MAP.txt")
    logger.info(f"OUTPUT WILL BE WRITTEN TO {output_file}")

    # Map BAM files to GFF
    n_bin = 1

    for bam_file in bam_file_list:
        bam_file_name = bam_file.split("/")[-1]

        # Map to stitched GFF
        mapped_out1 = os.path.join(mapped_folder, f"{stitched_gff_name}_{bam_file_name}_MAPPED.gff")
        cmd1 = [
            "rose2-bamToGFF",
            "-f", "1",
            "-e", "200",
            "-r",
            "-m", str(n_bin),
            "-b", bam_file,
            "-i", stitched_gff_file,
            "-o", mapped_out1,
        ]
        logger.info(f"Running: {' '.join(cmd1)}")
        subprocess.Popen(cmd1)

        # Map to original GFF
        mapped_out2 = os.path.join(mapped_folder, f"{input_name}_{bam_file_name}_MAPPED.gff")
        cmd2 = [
            "rose2-bamToGFF",
            "-f", "1",
            "-e", "200",
            "-r",
            "-m", str(n_bin),
            "-b", bam_file,
            "-i", input_gff_file,
            "-o", mapped_out2,
        ]
        logger.info(f"Running: {' '.join(cmd2)}")
        subprocess.Popen(cmd2)

    # Wait for mapping to complete
    logger.info("PAUSING TO MAP")
    time.sleep(10)

    # Check for mapping output
    output_done = False
    ticker = 0
    logger.info("WAITING FOR MAPPING TO COMPLETE. ELAPSED TIME (MIN):")

    while not output_done:
        output_done = True
        if ticker % 6 == 0:
            logger.info(f"{ticker * 5} minutes elapsed")
        ticker += 1

        if ticker == 144:  # 12 hours timeout
            logger.error("ERROR: OPERATION TIME OUT. MAPPING OUTPUT NOT DETECTED")
            return 1

        for bam_file in bam_file_list:
            bam_file_name = bam_file.split("/")[-1]
            mapped_out1 = os.path.join(
                mapped_folder, f"{stitched_gff_name}_{bam_file_name}_MAPPED.gff"
            )
            mapped_out2 = os.path.join(mapped_folder, f"{input_name}_{bam_file_name}_MAPPED.gff")

            for mapped_out in [mapped_out1, mapped_out2]:
                if not Path(mapped_out).exists():
                    output_done = False

        if output_done:
            break

        time.sleep(300)  # 5 minutes

    logger.info(f"MAPPING TOOK {ticker * 5} MINUTES")

    # Calculate density by region
    logger.info("BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS")
    map_collection(
        stitched_collection,
        reference_collection,
        bam_file_list,
        mapped_folder,
        output_file,
        ref_name=stitched_gff_name,
    )

    time.sleep(10)

    # Call R script for super-enhancer calling
    logger.info("CALLING AND PLOTTING SUPER-STITCHED PEAKS")

    if args.control:
        control_name = args.control.split("/")[-1]
    else:
        control_name = "NONE"

    r_script = Path(__file__).parent / "R" / "ROSE_callSuper.R"
    cmd = ["Rscript", str(r_script), out_folder, output_file, input_name, control_name]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # Call gene mapper
    time.sleep(60)
    super_table_file = os.path.join(out_folder, f"{input_name}_SuperStitched.table.txt")
    all_table_file = os.path.join(out_folder, f"{input_name}_AllStitched.table.txt")

    control_flag = ["-c"] if args.control else []

    if args.custom_genome:
        cmd1 = [
            "rose2-geneMapper",
            "--custom", args.custom_genome,
            "-i", super_table_file,
            "-r",  # FIXED: -r is a boolean flag, doesn't take a value
        ] + control_flag

        cmd2 = [
            "rose2-geneMapper",
            "--custom", args.custom_genome,
            "-i", all_table_file,
            "-r",  # FIXED: -r is a boolean flag, doesn't take a value
        ] + control_flag
    else:
        cmd1 = [
            "rose2-geneMapper",
            "-g", args.genome,
            "-i", super_table_file,
            "-r",  # FIXED: -r is a boolean flag, doesn't take a value
        ] + control_flag

        cmd2 = [
            "rose2-geneMapper",
            "-g", args.genome,
            "-i", all_table_file,
            "-r",  # FIXED: -r is a boolean flag, doesn't take a value
        ] + control_flag

    # Gene mapper for super-enhancers
    logger.info(f"Running: {' '.join(cmd1)}")
    subprocess.run(cmd1, check=True)

    # Gene mapper for all stitched peaks
    logger.info(f"Running: {' '.join(cmd2)}")
    subprocess.run(cmd2, check=True)

    logger.info("ROSE ANALYSIS COMPLETE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
