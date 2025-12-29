"""
ROSE2: Rank Ordering of Super-Enhancers

A computational pipeline for the identification and analysis of super-enhancers
from ChIP-seq data.

Original algorithm developed by the Young Lab at the Whitehead Institute.
Python 3 port by St. Jude Children's Research Hospital.
Modernization and PyPI packaging by Ming (Tommy) Tang.

This package provides tools to:
- Identify super-enhancers from ChIP-seq data
- Rank enhancer regions by signal density
- Map stitched enhancers to genes

References:
    Whyte et al. (2013) "Master Transcription Factors and Mediator Establish
    Super-Enhancers at Key Cell Identity Genes" Cell 153(2):307-319
    https://pmc.ncbi.nlm.nih.gov/articles/PMC3841062/
"""

__version__ = "2.0.0"
__author__ = "Richard Young Lab, St. Jude Children's Research Hospital, Ming Tang"

from rose2 import utils

__all__ = ["utils"]
