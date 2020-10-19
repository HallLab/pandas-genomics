try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata

from . import arrays, io, scalars
from .accessors import GenotypeAccessor

__version__ = importlib_metadata.version(__name__)

__all__ = [__version__, GenotypeAccessor, arrays, io, scalars]
