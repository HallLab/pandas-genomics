import numpy as np
import pandas as pd

from pandas_genomics.scalars import MISSING_IDX


class EncodingMixin:
    """
    Genotype Mixin containing functions for performing encoding
    """

    def encode_additive(self) -> pd.arrays.IntegerArray:
        """
        Additive Encoding

        - Number of copies of non-reference allele
        - np.nan when any alleles are missing

        Returns
        -------
        pd.arrays.IntegerArray
        """
        allele_sum = (self.allele_idxs != 0).sum(axis=1).astype("float")
        allele_sum[(self.allele_idxs == MISSING_IDX).any(axis=1)] = np.nan
        return allele_sum

    def encode_dominant(self) -> pd.arrays.IntegerArray:
        """
        Dominant Encoding

        - 0 for Homozygous Reference
        - 1 for any other case
        - np.nan when any alleles are missing

        Returns
        -------
        pd.arrays.IntegerArray
        """
        has_minor = (self.allele_idxs != 0).any(axis=1).astype("float")
        has_minor[(self.allele_idxs == MISSING_IDX).any(axis=1)] = np.nan
        return has_minor

    def encode_recessive(self) -> pd.arrays.IntegerArray:
        """
        Recessive Encoding

        - 1 for Homozygous Non-reference
        - 0 for anything else
        - np.nan when any alleles are missing

        Returns
        -------
        pd.arrays.IntegerArray
        """
        all_minor = (self.allele_idxs != 0).all(axis=1).astype("float")
        all_minor[(self.allele_idxs == MISSING_IDX).any(axis=1)] = np.nan
        return all_minor

    def encode_codominant(self) -> pd.arrays.Categorical:
        """
        This encodes the genotype into three categories.  When utilized in regression, this results in two variables
        due to dummy encoding- "Het" as 0 or 1 and "Hom" as 0 or 1.  0 in both indicates "Ref".

        - 'Ref' for Homozygous Reference
        - 'Het' for Heterozygous
        - 'Hom' for Homozygous Non-Reference
        - pd.NA for missing
        - Raises an error if ploidy is not 2

        Returns
        -------
        pd.arrays.Categorical
        """
        if self.dtype.variant.ploidy != 2:
            raise ValueError(
                "Codominant encoding can only be used with diploid genotypes"
            )

        allele_sum = (self.allele_idxs != 0).sum(axis=1)
        categories = ["Ref", "Het", "Hom"]
        result = pd.Categorical(
            values=[categories[n] if n in {0, 1, 2} else None for n in allele_sum],
            categories=categories,
            ordered=True,
        )
        result[(self.allele_idxs == MISSING_IDX).any(axis=1)] = None
        return result

    def encode_weighted(
        self,
        alpha_value: float,
        ref_allele: str,
        alt_allele: str,
        minor_allele_freq: float,
    ) -> pd.DataFrame:
        """Perform weighted (edge) encoding.

        Parameters
        ----------
        alpha_value: float
            Alpha value to use for heterozygous genotypes
        ref_allele: str
            Reference allele when the alpha was calculated
        alt_allele: str
            Alternate allele when the alpha was calculated
        minor_allele_freq: float
            Minor allele frequency when the alpha was calculated

        Returns
        -------
        pd.arrays.IntegerArray
            1 for Homozygous Alt
            alpha for Heterozygous Alt
            0 for Homozygous Ref
            np.NaN for any other genotype

        Notes
        -----
        See [1]_ for more information about weighted encoding.
        Encoding will be based on the provided ref and alt alleles, since the alpha value would be specific to them.
        In the future, if the existing minor allele frequency differs greatly from the provided value, a warning will be issued.

        References
        ----------
        .. [1] Hall, Molly A., et al.
               "Novel EDGE encoding method enhances ability to identify genetic interactions."
               PLoS genetics 17.6 (2021): e1009534.
        """
        # Get ref and alt allele idxes, throwing an error if the allele doesn't exist in the variant
        ref_allele_idx = self.variant.get_idx_from_allele(ref_allele)
        alt_allele_idx = self.variant.get_idx_from_allele(alt_allele)

        # TODO: Validate MAF.  Need to determine a reasonable warning threshold.

        encoded_values = pd.array(np.full(len(self), np.nan))
        encoded_values[(self.allele_idxs == ref_allele_idx).all(axis=1)] = 0.0
        encoded_values[(self.allele_idxs == alt_allele_idx).all(axis=1)] = 1.0
        encoded_values[
            (self.allele_idxs == sorted([ref_allele_idx, alt_allele_idx])).all(axis=1)
        ] = alpha_value

        return encoded_values
