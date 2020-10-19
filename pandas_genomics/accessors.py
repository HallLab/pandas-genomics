import pandas as pd

from pandas_genomics.arrays import GenotypeDtype


@pd.api.extensions.register_series_accessor("genotype")
class GenotypeAccessor:
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
    # In-place methods #
    ####################
    def set_reference(self, allele):
        self._array.set_reference(allele)

    ############
    # Encoding #
    ############
    def encode_additive(self):
        return pd.Series(
            data=self._array.encode_additive(),
            index=self._index,
            name=f"{self._array.variant.id}_{self._array.variant.alleles[1]}",
        )

    def encode_dominant(self):
        return pd.Series(
            data=self._array.encode_dominant(),
            index=self._index,
            name=f"{self._array.variant.id}_{self._array.variant.alleles[1]}",
        )

    def encode_recessive(self):
        return pd.Series(
            data=self._array.encode_recessive(),
            index=self._index,
            name=f"{self._array.variant.id}_{self._array.variant.alleles[1]}",
        )

    def encode_codominant(self):
        return pd.Series(
            data=self._array.encode_codominant(),
            index=self._index,
            name=f"{self._array.variant.id}_{self._array.variant.alleles[1]}",
        )
