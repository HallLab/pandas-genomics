import pandas as pd

from pandas_genomics.arrays import GenotypeDtype


@pd.api.extensions.register_dataframe_accessor("genomics")
class GenotypeDataframeAccessor:
    """
    DataFrame accessor for GenotypeArray methods
    """
    def __init__(self, pandas_obj):
        for colname in pandas_obj.columns:
            test = pandas_obj[colname]
            if not GenotypeDtype.is_dtype(pandas_obj[colname].values.dtype):
                raise AttributeError(
                    f"Incompatible datatype: column {colname}  is '{pandas_obj[colname].dtype}',"
                    f" but must be a GenotypeDtype"
                )
        self._obj = pandas_obj

    @property
    def variant_info(self):
        """Return a DataFrame with variant info indexed by the column name"""
        return pd.DataFrame.from_dict({
            colname: self._obj[colname].genomics.variant_info
            for colname in self._obj.columns
        }, orient="index")
