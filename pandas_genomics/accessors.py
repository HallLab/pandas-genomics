import pandas as pd

from pandas_genomics.arrays import GenotypeDtype


@pd.api.extensions.register_series_accessor("genotype")
class GenotypeAccessor:
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
    def variant_score(self):
        """Retrieve the variant score as a float

        Returns
        -------
        variant_score: float
        """
        if self._array.variant.score is None:
            return float("NaN")
        else:
            return float(self._array.variant.score)

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
