from collections import Counter
from typing import Optional, List

import numpy as np
import pandas as pd

from .utils import generate_weighted_encodings
from pandas_genomics.arrays import GenotypeDtype


@pd.api.extensions.register_dataframe_accessor("genomics")
class GenotypeDataframeAccessor:
    """
    DataFrame accessor for GenotypeArray methods

    .. code-block:: python

        df.genomics.variant_info
        df.genomics.encode_additive()
    """

    def __init__(self, pandas_obj):
        if not pandas_obj.dtypes.apply(lambda dt: GenotypeDtype.is_dtype(dt)).any():
            raise AttributeError(
                "Incompatible datatypes: at least one column must be a GenotypeDtype."
            )
        id_counts = Counter(
            [
                s.genomics.variant.id
                for _, s in pandas_obj.iteritems()
                if GenotypeDtype.is_dtype(s)
            ]
        )
        if len(id_counts) < len(pandas_obj.select_dtypes([GenotypeDtype]).columns):
            duplicates = [(k, v) for k, v in id_counts.items() if v >= 2]
            raise AttributeError(
                f"Duplicate Variant IDs.  Column names may differ from variant IDs, but variant IDs must be unique.\n\tDuplicates: "
                + ", ".join([f"{dupe} ({count:,})" for dupe, count in duplicates])
            )
        self._obj = pandas_obj

    ######################
    # Variant Properties #
    ######################
    # These methods generally only return a result for each GenotypeArray column, ignoring other columns

    @property
    def variant_info(self) -> pd.DataFrame:
        """Return a DataFrame with variant info indexed by the column name (one row per GenotypeArray)"""
        genotypes = self._obj.select_dtypes([GenotypeDtype])
        return pd.DataFrame.from_dict(
            {
                colname: series.genomics.variant_info
                for colname, series in genotypes.iteritems()
            },
            orient="index",
        )

    @property
    def maf(self):
        """Return the minor allele frequency of each variant

        See :py:attr:`GenotypeArray.maf`"""
        genotypes = self._obj.select_dtypes([GenotypeDtype])
        return genotypes.apply(lambda col: col.genomics.maf)

    @property
    def hwe_pval(self):
        """Return the probability that the samples are in HWE

        See :py:attr:`GenotypeArray.hwe_pval`"""
        genotypes = self._obj.select_dtypes([GenotypeDtype])
        return genotypes.apply(lambda col: col.genomics.hwe_pval)

    ############
    # Encoding #
    ############
    # These methods generally return encoded values for any GenotypeArray columns without modifying other columns

    def encode_additive(self) -> pd.DataFrame:
        """Additive encoding of genotypes.

        See :meth:`GenotypeArray.encode_additive`

        Returns
        -------
        pd.DataFrame
        """
        return pd.concat(
            [
                s.genomics.encode_additive() if GenotypeDtype.is_dtype(s) else s
                for _, s in self._obj.iteritems()
            ],
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
            [
                s.genomics.encode_dominant() if GenotypeDtype.is_dtype(s) else s
                for _, s in self._obj.iteritems()
            ],
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
            [
                s.genomics.encode_recessive() if GenotypeDtype.is_dtype(s) else s
                for _, s in self._obj.iteritems()
            ],
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
            [
                s.genomics.encode_codominant() if GenotypeDtype.is_dtype(s) else s
                for _, s in self._obj.iteritems()
            ],
            axis=1,
        )

    def encode_weighted(self, encoding_info: pd.DataFrame) -> pd.DataFrame:
        """Weighted (EDGE) encoding of genotypes.

        See :meth:`GenotypeArray.encode_weighted`

        Parameters
        ----------
        encoding_info: pd.DataFrame
            columns:
                - Variant ID - used to match variants
                - Alpha Value - used for heterozygous genotypes
                - Ref Allele - which allele is considered reference
                - Alt Allele - which allele is considered alternate
                - Minor Allele Frequency - MAF of data used during calculation of alpha values

        Returns
        -------
        pd.DataFrame
        """
        # Validate the input DataFrame
        for required_col in [
            "Variant ID",
            "Alpha Value",
            "Ref Allele",
            "Alt Allele",
            "Minor Allele Frequency",
        ]:
            if required_col not in list(encoding_info):
                raise ValueError(
                    f"Missing one or more required columns in the encoding info: `{required_col}`"
                )
        id_counts = encoding_info["Variant ID"].value_counts()
        if sum(id_counts > 1):
            raise ValueError(
                f"Duplicate IDs: {', '.join([v for v in id_counts[id_counts>1].index])}"
            )

        # Rename the columns to match parameter names for simplicity
        encoding_info = encoding_info.rename(
            columns={
                "Alpha Value": "alpha_value",
                "Ref Allele": "ref_allele",
                "Alt Allele": "alt_allele",
                "Minor Allele Frequency": "minor_allele_freq",
            }
        )

        # Convert the encoding info into a Dict("Variant ID" = {param names : param values})
        encoding_info = {
            d["Variant ID"]: {k: v for k, v in d.items() if k != "Variant ID"}
            for d in encoding_info.to_dict(orient="records")
        }

        # Log messages for any warnings
        warnings = dict()

        # Process each variant
        results = []
        for _, s in self._obj.iteritems():
            if not GenotypeDtype.is_dtype(s):
                results.append(s)
                continue
            info = encoding_info.get(s.array.variant.id, None)
            if info is None:
                warnings[
                    s.array.variant.id
                ] = "No matching information found in the encoding data"
                continue
            elif (s.genomics.maf / info["minor_allele_freq"]) > 10e30:
                # TODO: replace this with a reasonable comparison to the data MAF.  For now it is an always-pass criteria
                warnings[
                    s.array.variant.id
                ] = f"Large MAF Difference: {s.genomics.maf} in sample, {info['minor_allele_freq']} in encoding data"
                continue
            else:
                try:
                    results.append(s.genomics.encode_weighted(**info))
                except Exception as e:
                    warnings[s.array.variant.id] = str(e)
        # Print Warnings
        if len(warnings) > 0:
            print(f"{len(warnings):,} Variables failed encoding")
            for var, warning in warnings.items():
                print(f"\t{var}: {warning}")
        # Concatenate results
        return pd.concat(results, axis=1)

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
            genotypes=self._obj.select_dtypes([GenotypeDtype]),
            data=data,
            outcome_variable=outcome_variable,
            covariates=covariates,
        )

    ###########
    # Filters #
    ###########
    # These methods drop genotypes that fail the filter, ignoring other columns

    def filter_variants_maf(self, keep_min_freq: float = 0.01) -> pd.DataFrame:
        """
        Drop variants with a MAF less than the specified value (0.01 by default)
        """
        genotypes = self._obj.select_dtypes([GenotypeDtype])
        removed = genotypes.loc[:, genotypes.genomics.maf < keep_min_freq].columns
        return self._obj.drop(columns=removed)

    def filter_variants_hwe(self, cutoff: float = 0.05) -> pd.DataFrame:
        """
        Drop variants with a probability of HWE less than the specified value (0.05 by default).
        Keep np.nan results, which occur for non-diploid variants and insufficient sample sizes
        """
        genotypes = self._obj.select_dtypes([GenotypeDtype])
        genotype_hwe_pval = genotypes.genomics.hwe_pval
        removed = genotypes.loc[
            :, (genotype_hwe_pval < cutoff) & ~np.isnan(genotype_hwe_pval)
        ].columns
        return self._obj.drop(columns=removed)
