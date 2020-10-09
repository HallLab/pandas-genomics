"""
Input/Output
------------
The `io` module provides functions for loading and saving GenotypeArrays to common variant formats

.. autosummary::
     :toctree: io

     from_plink
"""

from .plink import from_plink

__all__ = ['from_plink', ]
