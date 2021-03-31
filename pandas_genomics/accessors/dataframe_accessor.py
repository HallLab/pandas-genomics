import pandas as pd

from pandas_genomics.arrays import GenotypeDtype


@pd.api.extensions.register_dataframe_accessor("genomics")
class GenotypeDataframeAccessor:
    """
    DataFrame accessor for GenotypeArray methods
    """

    def __init__(self, pandas_obj):
        for colname in pandas_obj.columns:
            if not GenotypeDtype.is_dtype(pandas_obj[colname].values.dtype):
                raise AttributeError(
                    f"Incompatible datatype: column {colname}  is '{pandas_obj[colname].dtype}',"
                    f" but must be a GenotypeDtype"
                )
        self._obj = pandas_obj

    ######################
    # Variant Properties #
    ######################
    @property
    def variant_info(self) -> pd.DataFrame:
        """Return a DataFrame with variant info indexed by the column name"""
        return pd.DataFrame.from_dict(
            {
                colname: self._obj[colname].genomics.variant_info
                for colname in self._obj.columns
            },
            orient="index",
        )

    ############
    # Encoding #
    ############
    def encode_additive(self) -> pd.DataFrame:
        """Additive encoding of genotypes.

        See :meth:`GenotypeArray.encode_additive`

        Returns
        -------
        pd.DataFrame
        """
        return pd.concat(
            [self._obj[col].genomics.encode_additive() for col in self._obj.columns],
            axis=1,
        )

    def encode_dominant(self) -> pd.DataFrame:
        """Dominant encoding of genotypes.

        See :meth:`GenotypeArray.encode_dominant`

        Returns
        -------
        pd.DataFrame
        """
        return pd.concat(
            [self._obj[col].genomics.encode_dominant() for col in self._obj.columns],
            axis=1,
        )

    def encode_recessive(self) -> pd.DataFrame:
        """Recessive encoding of genotypes.

        See :meth:`GenotypeArray.encode_recessive`

        Returns
        -------
        pd.DataFrame
        """
        return pd.concat(
            [self._obj[col].genomics.encode_recessive() for col in self._obj.columns],
            axis=1,
        )

    def encode_codominant(self) -> pd.DataFrame:
        """Codominant encoding of genotypes.

        See :meth:`GenotypeArray.encode_codominant`

        Returns
        -------
        pd.DataFrame
        """
        return pd.concat(
            [self._obj[col].genomics.encode_codominant() for col in self._obj.columns],
            axis=1,
        )
