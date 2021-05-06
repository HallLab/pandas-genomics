from typing import Optional

import numpy as np
import pandas as pd

from pandas_genomics.arrays import GenotypeDtype


@pd.api.extensions.register_dataframe_accessor("genomics")
class GenotypeDataframeAccessor:
    """
    DataFrame accessor for GenotypeArray methods
    """

    def __init__(self, pandas_obj):
        if not pandas_obj.dtypes.apply(lambda dt: GenotypeDtype.is_dtype(dt)).all():
            incorrect = pandas_obj.dtypes[
                ~pandas_obj.dtypes.apply(lambda dt: GenotypeDtype.is_dtype(dt))
            ]
            raise AttributeError(
                f"Incompatible datatypes: all columns must be a GenotypeDtype: {incorrect}"
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
                colname: series.genomics.variant_info
                for colname, series in self._obj.iteritems()
            },
            orient="index",
        )

    #########################
    # Calculated Properties #
    #########################
    @property
    def maf(self):
        """Return the minor allele frequency

        See :py:attr:`GenotypeArray.maf`"""
        return self._obj.apply(lambda col: col.genomics.maf)

    @property
    def hwe_pval(self):
        """Return the probability that the samples are in HWE

        See :py:attr:`GenotypeArray.hwe_pval`"""
        return self._obj.apply(lambda col: col.genomics.hwe_pval)

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
            [s.genomics.encode_additive() for _, s in self._obj.iteritems()], axis=1
        )

    def encode_dominant(self) -> pd.DataFrame:
        """Dominant encoding of genotypes.

        See :meth:`GenotypeArray.encode_dominant`

        Returns
        -------
        pd.DataFrame
        """
        return pd.concat(
            [s.genomics.encode_dominant() for _, s in self._obj.iteritems()], axis=1
        )

    def encode_recessive(self) -> pd.DataFrame:
        """Recessive encoding of genotypes.

        See :meth:`GenotypeArray.encode_recessive`

        Returns
        -------
        pd.DataFrame
        """
        return pd.concat(
            [s.genomics.encode_recessive() for _, s in self._obj.iteritems()], axis=1
        )

    def encode_codominant(self) -> pd.DataFrame:
        """Codominant encoding of genotypes.

        See :meth:`GenotypeArray.encode_codominant`

        Returns
        -------
        pd.DataFrame
        """
        return pd.concat(
            [s.genomics.encode_codominant() for _, s in self._obj.iteritems()], axis=1
        )

    ###########
    # Filters #
    ###########
    def filter_variants_maf(self, keep_min_freq: float = 0.01) -> pd.DataFrame:
        """
        Drop variants with a MAF less than the specified value (0.01 by default)
        """
        return self._obj.loc[:, self._obj.genomics.maf >= keep_min_freq]

    def filter_variants_hwe(self, cutoff: float = 0.05) -> pd.DataFrame:
        """
        Drop variants with a probability of HWE less than the specified value (0.05 by default).
        Keep np.nan results, which occur for non-diploid variants and insufficient sample sizes
        """
        return self._obj.loc[
            :,
            (self._obj.genomics.hwe_pval >= cutoff)
            | (np.isnan(self._obj.genomics.hwe_pval)),
        ]
