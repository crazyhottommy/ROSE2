#!/usr/bin/env python3
"""
ROSE2 Command Line Interface

Main entry point for the ROSE2 package providing access to all sub-commands.

Original: Young Lab, Whitehead Institute
Python 3 Port: St. Jude Children's Research Hospital
Modernization: Ming Tang
"""

import argparse
import sys
from typing import List, Optional

from rose2 import __version__


def main(argv: Optional[List[str]] = None) -> int:
    """Main CLI entry point.

    Args:
        argv: Command line arguments (for testing)

    Returns:
        Exit code
    """
    parser = argparse.ArgumentParser(
        prog="rose2",
        description="ROSE2: Rank Ordering of Super-Enhancers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full ROSE analysis
  rose2 main -i peaks.gff -r sample.bam -g HG19 -o output/

  # Map BAM to GFF
  rose2 bamToGFF -b sample.bam -i regions.gff -o output.gff

  # Map enhancers to genes
  rose2 geneMapper -i enhancers.txt -g HG19 -o output/

For more information, see: https://github.com/stjude/ROSE2
        """,
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    subparsers = parser.add_subparsers(
        title="commands",
        description="Available ROSE2 commands",
        dest="command",
        help="Use 'rose2 <command> --help' for command-specific help",
    )

    # Main command
    subparsers.add_parser(
        "main",
        help="Run full ROSE analysis pipeline",
        add_help=False,
    )

    # BamToGFF command
    subparsers.add_parser(
        "bamToGFF",
        help="Map BAM reads to GFF regions",
        add_help=False,
    )

    # GeneMapper command
    subparsers.add_parser(
        "geneMapper",
        help="Map enhancers to genes",
        add_help=False,
    )

    args, remaining = parser.parse_known_args(argv)

    if not args.command:
        parser.print_help()
        return 0

    # Import and call the appropriate module
    if args.command == "main":
        from rose2.rose_main import main as rose_main
        return rose_main(remaining)

    elif args.command == "bamToGFF":
        from rose2.rose_bamToGFF import main as bam_to_gff_main
        return bam_to_gff_main(remaining)

    elif args.command == "geneMapper":
        from rose2.rose_geneMapper import main as gene_mapper_main
        return gene_mapper_main(remaining)

    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
