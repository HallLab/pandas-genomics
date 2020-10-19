from . import arrays, io, scalars
from .accessors import GenotypeAccessor

# Simple version tracking for now until Poetry has a solution
__version__ = "v0.3.0"

__all__ = [__version__, GenotypeAccessor, arrays, io, scalars]
