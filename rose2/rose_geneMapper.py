#!/usr/bin/env python3
"""
ROSE Gene Mapper

Maps enhancer regions to genes based on proximity and overlap.
Creates two outputs: a gene-mapped region table where each row is an enhancer,
and a gene table where each row is a gene.

Original: Young Lab, Whitehead Institute
Python 3 Port: St. Jude Children's Research Hospital
Modernization: Ming Tang
"""

import argparse
import logging
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, List, Optional, Tuple

from rose2 import utils

logger = logging.getLogger(__name__)


def map_enhancer_to_gene(
    annot_file: str,
    enhancer_file: str,
    transcribed_file: str = "",
    unique_genes: bool = True,
    by_refseq: bool = False,
    subtract_input: bool = False,
) -> Tuple[List[List[Any]], List[List[Any]], List[List[Any]]]:
    """Map enhancers to genes based on proximity and overlap.

    Args:
        annot_file: Path to annotation file
        enhancer_file: Path to enhancer table file
        transcribed_file: Path to file with transcribed genes (optional)
        unique_genes: If True, reduce to gene name only
        by_refseq: If True, output by RefSeq ID instead of gene name
        subtract_input: If True, subtract input signal from sample signal

    Returns:
        Tuple of (enhancer-to-gene table, gene-to-enhancer table, signal-with-genes table)
    """
    logger.info("Creating start dictionary")
    start_dict = utils.make_start_dict(annot_file)

    logger.info("Parsing enhancer table")
    enhancer_table = utils.parse_table(enhancer_file, "\t")

    # Get list of transcribed genes
    if len(transcribed_file) > 0:
        transcribed_table = utils.parse_table(transcribed_file, "\t")
        transcribed_genes = [line[1] for line in transcribed_table]
    else:
        transcribed_genes = list(start_dict.keys())

    # Create transcript collection
    logger.info("Making transcript collection")
    transcribed_collection = utils.make_transcript_collection(
        annot_file, 0, 0, 500, transcribed_genes
    )

    # Create TSS collection
    logger.info("Making TSS collection")
    tss_loci = []
    for gene_id in transcribed_genes:
        tss_loci.append(utils.make_tss_locus(gene_id, start_dict, 0, 0))

    tss_collection = utils.LocusCollection(tss_loci, 50)

    # Initialize gene dictionaries
    gene_dict: Dict[str, DefaultDict[str, List[str]]] = {
        "overlapping": defaultdict(list),
        "proximal": defaultdict(list),
        "enhancerString": defaultdict(list),
    }
    overall_gene_list: List[str] = []

    # Set up output tables
    enhancer_to_gene_table = [
        enhancer_table[5][0:6] + ["OVERLAP_GENES", "PROXIMAL_GENES", "CLOSEST_GENE"]
        + enhancer_table[5][-2:]
    ]

    gene_to_enhancer_table = [["GENE_NAME", "REFSEQ_ID", "PROXIMAL_STITCHED_PEAKS"]]

    signal_with_genes = [["GENE_NAME", "REFSEQ_ID", "PROXIMAL_STITCHED_PEAKS", "SIGNAL"]]

    # Process each enhancer
    for line in enhancer_table[6:]:
        enhancer_string = f"{line[1]}:{line[2]}-{line[3]}"
        enhancer_signal = int(float(line[6]))
        if subtract_input:
            enhancer_signal = int(float(line[6]) - float(line[7]))

        enhancer_locus = utils.Locus(line[1], line[2], line[3], ".", line[0])

        # Find overlapping genes (genes whose transcript overlaps the enhancer)
        overlapping_loci = transcribed_collection.get_overlap(enhancer_locus, "both")
        overlapping_genes = [locus.ID() for locus in overlapping_loci]

        # Find proximal genes (TSS within 50kb of enhancer boundary)
        proximal_loci = tss_collection.get_overlap(
            utils.make_search_locus(enhancer_locus, 50000, 50000), "both"
        )
        proximal_genes = [locus.ID() for locus in proximal_loci]

        # OPTIMIZED: Instead of searching 50Mb window (2M window checks = 6 hours!),
        # use smarter logic:
        # - If we have overlapping/proximal genes (≤50kb), they will ALWAYS be closest
        # - Only search further if NO nearby genes found
        # - Use graduated search to avoid checking millions of windows
        distal_genes = []

        if not overlapping_genes and not proximal_genes:
            # No nearby genes - need to search further for closest gene
            # Use graduated distances to minimize window checks
            search_distances = [500000, 2000000, 10000000]  # 500kb, 2Mb, 10Mb

            for dist in search_distances:
                distal_loci = tss_collection.get_overlap(
                    utils.make_search_locus(enhancer_locus, dist, dist), "both"
                )
                distal_genes = [locus.ID() for locus in distal_loci]
                if distal_genes:
                    break  # Found genes, stop expanding
            # If still no genes after 10Mb, closest_gene will be empty (handled below)

        # Uniquify gene lists
        overlapping_genes = utils.uniquify(overlapping_genes)
        proximal_genes = utils.uniquify(proximal_genes)
        distal_genes = utils.uniquify(distal_genes)
        all_enhancer_genes = overlapping_genes + proximal_genes + distal_genes

        # Remove overlapping genes from proximal list
        for ref_id in overlapping_genes:
            if ref_id in proximal_genes:
                proximal_genes.remove(ref_id)

        # Remove proximal genes from distal list
        for ref_id in proximal_genes:
            if ref_id in distal_genes:
                distal_genes.remove(ref_id)

        # Find closest gene
        if len(all_enhancer_genes) == 0:
            closest_gene = ""
        else:
            # Calculate enhancer center
            enhancer_center = (int(line[2]) + int(line[3])) / 2

            # Get absolute distance to enhancer center
            dist_list = [
                abs(enhancer_center - start_dict[gene_id]["start"][0])
                for gene_id in all_enhancer_genes
            ]
            closest_gene_id = all_enhancer_genes[dist_list.index(min(dist_list))]

            if by_refseq:
                closest_gene = closest_gene_id
            else:
                closest_gene = start_dict[closest_gene_id]["name"]

        # Write enhancer table row
        new_enhancer_line = line[0:6]

        if by_refseq:
            new_enhancer_line.append(",".join(utils.uniquify(overlapping_genes)))
            new_enhancer_line.append(",".join(utils.uniquify(proximal_genes)))
            new_enhancer_line.append(closest_gene)
        else:
            new_enhancer_line.append(
                ",".join(utils.uniquify([start_dict[x]["name"] for x in overlapping_genes]))
            )
            new_enhancer_line.append(
                ",".join(utils.uniquify([start_dict[x]["name"] for x in proximal_genes]))
            )
            new_enhancer_line.append(closest_gene)

        # Write signal with genes table
        if closest_gene:
            if by_refseq:
                signal_with_genes.append(
                    [start_dict[closest_gene]["name"], closest_gene, enhancer_string, enhancer_signal]
                )
            else:
                # Find the RefSeq ID for the gene name
                closest_refseq = ""
                for gene_id in all_enhancer_genes:
                    if start_dict[gene_id]["name"] == closest_gene:
                        closest_refseq = gene_id
                        break
                signal_with_genes.append([closest_gene, closest_refseq, enhancer_string, enhancer_signal])

        new_enhancer_line += line[-2:]
        enhancer_to_gene_table.append(new_enhancer_line)

        # Add genes to overall list and dictionaries
        overall_gene_list += overlapping_genes
        for ref_id in overlapping_genes:
            gene_dict["overlapping"][ref_id].append(enhancer_string)

        overall_gene_list += proximal_genes
        for ref_id in proximal_genes:
            gene_dict["proximal"][ref_id].append(enhancer_string)

    # Create gene-to-enhancer table
    overall_gene_list = utils.uniquify(overall_gene_list)

    name_order = utils.order([start_dict[x]["name"] for x in overall_gene_list])
    used_names: List[str] = []

    for i in name_order:
        ref_id = overall_gene_list[i]
        gene_name = start_dict[ref_id]["name"]

        if gene_name in used_names and unique_genes:
            continue
        else:
            used_names.append(gene_name)

        prox_enhancers = gene_dict["proximal"][ref_id] + gene_dict["overlapping"][ref_id]
        new_line = [gene_name, ref_id, ",".join(prox_enhancers)]
        gene_to_enhancer_table.append(new_line)

    # Re-sort enhancer-to-gene table by signal
    enhancer_order = utils.order([int(line[-2]) for line in enhancer_to_gene_table[1:]])
    sorted_table = [enhancer_to_gene_table[0]]
    for i in enhancer_order:
        sorted_table.append(enhancer_to_gene_table[i + 1])

    return sorted_table, gene_to_enhancer_table, signal_with_genes


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
        description="Map ROSE enhancers to genes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required arguments
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input", required=True, help="ROSE ranked enhancer or super-enhancer file"
    )

    # Genome arguments (one required)
    genome_group = parser.add_mutually_exclusive_group(required=True)
    genome_group.add_argument(
        "-g",
        "--genome",
        choices=["MM9", "MM8", "MM10", "HG18", "HG19", "HG38"],
        help="Genome build",
    )
    genome_group.add_argument(
        "--custom", dest="custom_genome", help="Custom genome annotation .ucsc file"
    )

    # Optional arguments
    parser.add_argument("-l", "--list", dest="gene_list", help="Gene list to filter through")
    parser.add_argument(
        "-o",
        "--out",
        help="Output folder (default: same folder as input file)",
    )
    parser.add_argument(
        "-r",
        "--refseq",
        action="store_true",
        default=False,
        help="Output by RefSeq ID instead of common name",
    )
    parser.add_argument(
        "-c",
        "--control",
        action="store_true",
        default=False,
        help="Subtract input from sample signal",
    )

    args = parser.parse_args(argv)

    # Validate input file
    if not Path(args.input).exists():
        logger.error(f"Input file not found: {args.input}")
        return 1

    enhancer_file = args.input

    # Set output folder
    if args.out:
        out_folder = utils.format_folder(args.out, True)
    else:
        out_folder = str(Path(enhancer_file).parent) + "/"

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

    # Get transcribed gene list
    transcribed_file = args.gene_list if args.gene_list else ""

    # Map enhancers to genes
    enhancer_to_gene_table, gene_to_enhancer_table, with_genes_table = map_enhancer_to_gene(
        annot_file,
        enhancer_file,
        unique_genes=True,
        by_refseq=args.refseq,
        subtract_input=args.control,
        transcribed_file=transcribed_file,
    )

    # Get enhancer filename
    enhancer_file_name = Path(enhancer_file).stem

    # Write output files
    out1 = os.path.join(out_folder, f"{enhancer_file_name}_REGION_TO_GENE.txt")
    utils.unparse_table(enhancer_to_gene_table, out1, "\t")
    logger.info(f"Wrote enhancer-to-gene table: {out1}")

    out2 = os.path.join(out_folder, f"{enhancer_file_name}_GENE_TO_REGION.txt")
    utils.unparse_table(gene_to_enhancer_table, out2, "\t")
    logger.info(f"Wrote gene-to-enhancer table: {out2}")

    out3 = os.path.join(out_folder, f"{enhancer_file_name}.table_withGENES.txt")
    utils.unparse_table(with_genes_table, out3, "\t")
    logger.info(f"Wrote signal with genes table: {out3}")

    logger.info("Gene mapping complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
