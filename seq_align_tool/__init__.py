"""
生物序列比对工具核心包
"""

from .fasta_reader import read_fasta, read_fasta_directory
from .scoring import ScoringMatrix
from .alignment import (
    AlignmentResult,
    needleman_wunsch,
    smith_waterman
)
from .visualization import (
    color_alignment,
    print_alignment
)
from .export import export_alignment
from .batch import batch_alignment, export_batch_results
from .conservation import conservation_analysis
from .genbank import fetch_genbank

__all__ = [
    'read_fasta',
    'read_fasta_directory',
    'ScoringMatrix',
    'AlignmentResult',
    'needleman_wunsch',
    'smith_waterman',
    'color_alignment',
    'print_alignment',
    'export_alignment',
    'batch_alignment',
    'export_batch_results',
    'conservation_analysis',
    'fetch_genbank'
]

__version__ = '1.0.0'
