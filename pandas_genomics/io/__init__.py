"""
``pandas_genomics.io``
======================
This module provides functions for loading and saving GenotypeArrays to common variant formats

Loading data
------------
   from_plink

"""

from .plink import from_plink

__all__ = ['from_plink', ]
