from .dtypes.scalars import Genotype, Variant
from .dtypes.genotype import GenotypeDtype, GenotypeArray

__all__ = [
    'Genotype',
    'Variant',
    'GenotypeDtype',
    'GenotypeArray'
]

# Simple version tracking for now until Poetry has a solution
__version__ = "v0.1.0"