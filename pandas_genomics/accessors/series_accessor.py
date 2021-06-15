import pandas as pd

from pandas_genomics.arrays import GenotypeDtype


@pd.api.extensions.register_series_accessor("genomics")
class GenotypeSeriesAccessor:
    """
    Series accessor for GenotypeArray methods
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
            data=self._array.encode_additive(),
            index=self._index,
            name=f"{self._array.variant.id}_{self._array.variant.alleles[1]}",
        )

    def encode_dominant(self) -> pd.Series:
        """Dominant encoding of genotypes.

        See :meth:`GenotypeArray.encode_dominant`

        Returns
        -------
        pd.Series
        """
        return pd.Series(
            data=self._array.encode_dominant(),
            index=self._index,
            name=f"{self._array.variant.id}_{self._array.variant.alleles[1]}",
        )

    def encode_recessive(self) -> pd.Series:
        """Recessive encoding of genotypes.

        See :meth:`GenotypeArray.encode_recessive`

        Returns
        -------
        pd.Series
        """
        return pd.Series(
            data=self._array.encode_recessive(),
            index=self._index,
            name=f"{self._array.variant.id}_{self._array.variant.alleles[1]}",
        )

    def encode_codominant(self) -> pd.Series:
        """Codominant encoding of genotypes.

        See :meth:`GenotypeArray.encode_codominant`

        Returns
        -------
        pd.Series
        """
        return pd.Series(
            data=self._array.encode_codominant(),
            index=self._index,
            name=f"{self._array.variant.id}_{self._array.variant.alleles[1]}",
        )

    def encode_weighted(self,
                        alpha_value: float,
                        ref_allele: str,
                        alt_allele: str,
                        minor_allele_freq: float) -> pd.Series:
        """Weighted (edge) encoding of genotypes.

        See :meth:`GenotypeArray.encode_weighted`

        Returns
        -------
        pd.Series
        """
        return pd.Series(
            data=self._array.encode_weighted(alpha_value, ref_allele, alt_allele, minor_allele_freq),
            index=self._index,
            name=f"{self._array.variant.id}_{alt_allele}",
        )

    ##############
    # QC Methods #
    ##############
