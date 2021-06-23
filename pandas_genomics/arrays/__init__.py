"""
.. currentmodule:: pandas_genomics.arrays

Arrays
------

This module contains ExtensionArrays and their corresponding ExtensionDtypes

  .. autosummary::
     :toctree: arrays

     GenotypeDtype
     GenotypeArray

Specialized methods are added to the GenotypeArray using Mixins:

  .. autosummary::
     : toctree: arrays

     encoding_mixin.EncodingMixin
     info_mixin.InfoMixin

"""

from .genotype_array import GenotypeDtype, GenotypeArray

__all__ = ["GenotypeDtype", "GenotypeArray"]
