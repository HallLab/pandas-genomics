from typing import Optional, List

import pandas as pd

from .utils import generate_weighted_encodings
from pandas_genomics.arrays import GenotypeDtype


@pd.api.extensions.register_series_accessor("genomics")
class GenotypeSeriesAccessor:
    """
    Series accessor for GenotypeArray methods


    .. code-block:: python

        s.genomics.variant_info
        s.genomics.encode_additive()
    """

    def __init__(self, obj):
        if not GenotypeDtype.is_dtype(obj.dtype):
            raise AttributeError(
                f"Incompatible datatype ({obj.dtype}), must be a GenotypeDtype"
            )
        self._array = obj.array
        self._index = obj.index
        self._name = obj.name

    def _wrap_method(self, method, *args, **kwargs):
        return pd.Series(method(*args, **kwargs), self.index, name=self.name)

    ####################
    # Variant Properties
    ####################
    @property
    def variant(self):
        """Retrieve the variant object

        Returns
        -------
        variant: Variant
        """
        return self._array.variant

    @property
    def variant_info(self):
        """Retrieve the variant as a pandas Series

        Returns
        -------
        variant: pd.Series"""
        return pd.Series(self._array.variant.as_dict(), name=self._name)

    #######################
    # Genotype Properties #
    #######################
    @property
    def gt_scores(self):
        """Return an array of genotype scores as float values
        np.nan when the score is missing
        """
        return self._array.gt_scores

    #########################
    # Calculated Properties #
    #########################
    @property
    def maf(self):
        """Return the minor allele frequency

        See :py:attr:`GenotypeArray.maf`"""
        return self._array.maf

    #########################
    # Calculated Properties #
    #########################
    @property
    def hwe_pval(self):
        """Return the probability that the samples are in HWE

        See :py:attr:`GenotypeArray.hwe_pval`"""
        return self._array.hwe_pval

    ####################
    # In-place methods #
    ####################
    def set_reference(self, allele) -> None:
        """Change the allele reference variant.

        See :meth:`GenotypeArray.set_reference`

        Parameters
        ----------
        allele: str
            Must match an allele already in the variant

        Returns
        -------
        None
        """
        self._array.set_reference(allele)

    ############
    # Encoding #
    ############
    def encode_additive(self) -> pd.Series:
        """Additive encoding of genotypes.

        See :meth:`GenotypeArray.encode_additive`

        Returns
        -------
        pd.Series
        """
        return pd.Series(
            data=self._array.encode_additive(), index=self._index, name=self._name
        )

    def encode_dominant(self) -> pd.Series:
        """Dominant encoding of genotypes.

        See :meth:`GenotypeArray.encode_dominant`

        Returns
        -------
        pd.Series
        """
        return pd.Series(
            data=self._array.encode_dominant(), index=self._index, name=self._name
        )

    def encode_recessive(self) -> pd.Series:
        """Recessive encoding of genotypes.

        See :meth:`GenotypeArray.encode_recessive`

        Returns
        -------
        pd.Series
        """
        return pd.Series(
            data=self._array.encode_recessive(), index=self._index, name=self._name
        )

    def encode_codominant(self) -> pd.Series:
        """Codominant encoding of genotypes.

        See :meth:`GenotypeArray.encode_codominant`

        Returns
        -------
        pd.Series
        """
        return pd.Series(
            data=self._array.encode_codominant(), index=self._index, name=self._name
        )

    def encode_weighted(
        self,
        alpha_value: float,
        ref_allele: str,
        alt_allele: str,
        minor_allele_freq: float,
    ) -> pd.Series:
        """Weighted (edge) encoding of genotypes.

        See :meth:`GenotypeArray.encode_weighted`

        Returns
        -------
        pd.Series
        """
        return pd.Series(
            data=self._array.encode_weighted(
                alpha_value, ref_allele, alt_allele, minor_allele_freq
            ),
            index=self._index,
            name=self._name,
        )

    def generate_weighted_encodings(
        self,
        data: pd.DataFrame,
        outcome_variable: str,
        covariates: Optional[List[str]] = None,
    ):
        """
        Calculate alpha values to be used in weighted encoding

        Parameters
        ----------
        data:
            Data to be used in the regression, including the outcome and covariates
        outcome_variable:
            The variable to be used as the output (y) of the regression
        covariates:
            Other variables to be included in the regression formula

        Returns
        -------
        Dict
          Variant ID: str
          Alpha Value - used for heterozygous genotypes
          Ref Allele - which allele is considered reference
          Alt Allele - which allele is considered alternate
          Minor Allele Frequency - MAF of data used during calculation of alpha values

        Notes
        -----
        See [1]_ for more information about weighted encoding.

        References
        ----------
        .. [1] Hall, Molly A., et al.
               "Novel EDGE encoding method enhances ability to identify genetic interactions."
               PLoS genetics 17.6 (2021): e1009534.
        """
        return generate_weighted_encodings(
            genotypes=pd.Series(self._array, name=self._name, index=self._index),
            data=data,
            outcome_variable=outcome_variable,
            covariates=covariates,
        )

    ##############
    # QC Methods #
    ##############
