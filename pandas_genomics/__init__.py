from .scalars import Genotype, Variant
from .arrays.genotype_array import GenotypeDtype, GenotypeArray

__all__ = [
    'Genotype',
    'Variant',
    'GenotypeDtype',
    'GenotypeArray'
]

# Simple version tracking for now until Poetry has a solution
__version__ = "v0.1.0"
