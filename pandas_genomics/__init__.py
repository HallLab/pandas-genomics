from .scalars import Genotype, Variant
from .arrays import GenotypeDtype, GenotypeArray
import io

__all__ = [
    'Genotype',
    'Variant',
    'GenotypeDtype',
    'GenotypeArray',
    'io'
]

# Simple version tracking for now until Poetry has a solution
__version__ = "v0.2.0"
