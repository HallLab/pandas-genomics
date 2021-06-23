"""
Accessors
---------

This module contains 'genomics' accessors for DataFrames and Series

  .. autosummary::
     :toctree: accessors

     GenotypeSeriesAccessor
     GenotypeDataframeAccessor

"""

from .series_accessor import GenotypeSeriesAccessor
from .dataframe_accessor import GenotypeDataframeAccessor

__all__ = ["GenotypeSeriesAccessor", "GenotypeDataframeAccessor"]
