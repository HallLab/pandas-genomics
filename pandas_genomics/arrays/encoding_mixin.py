import numpy as np
import pandas as pd

from pandas_genomics.scalars import MISSING_IDX


class EncodingMixin:
    def encode_additive(self) -> pd.arrays.IntegerArray:
        """
        Returns
        -------
        pd.arrays.IntegerArray
            Number of copies of the minor allele
            pd.NA when any alleles are missing
            Raises ValueError if there is more than 1 alternate allele
        """
        # TODO: Return multiple arrays for multiple alternate alleles?
        if len(self.variant.alleles) > 2:
            raise ValueError("Additive encoding can only be used with one allele")

        allele_sum = self.allele_idxs.sum(axis=1).astype("float")
        allele_sum[(self.allele_idxs == MISSING_IDX).any(axis=1)] = np.nan
        result = pd.array(data=allele_sum, dtype="UInt8")
        return result

    def encode_dominant(self) -> pd.arrays.IntegerArray:
        """
        Returns
        -------
        pd.arrays.IntegerArray
            0 for no copies of the minor allele
            1 for any copies of the minor allele
            pd.NA when any alleles are missing
            Raises an error if there is more than 1 alternate allele
        """
        # TODO: Return multiple arrays for multiple alternate alleles?
        if len(self.variant.alleles) > 2:
            raise ValueError("Dominant encoding can only be used with one allele")

        has_minor = (self.allele_idxs == 1).any(axis=1).astype("float")
        has_minor[(self.allele_idxs == MISSING_IDX).any(axis=1)] = np.nan
        result = pd.array(data=has_minor, dtype="UInt8")
        return result

    def encode_recessive(self) -> pd.arrays.IntegerArray:
        """
        Returns
        -------
        pd.arrays.IntegerArray
            1 for Homozygous Alt
            0 for anything else
            pd.NA when any alleles are missing
            Raises an error if there is more than 1 alternate allele
        """
        # TODO: Return multiple arrays for multiple alternate alleles?
        if len(self.variant.alleles) > 2:
            raise ValueError("Recessive encoding can only be used with one allele")

        all_minor = (self.allele_idxs == 1).all(axis=1).astype("float")
        all_minor[(self.allele_idxs == MISSING_IDX).any(axis=1)] = np.nan
        result = pd.array(data=all_minor, dtype="UInt8")
        return result

    def encode_codominant(self) -> pd.arrays.Categorical:
        """
        This encodes the genotype into three categories.  When utilized in regression, this results in two variables
        due to dummy encoding- "Het" as 0 or 1 and "Hom" as 0 or 1.  0 in both indicates "Ref".

        Returns
        -------
        pd.arrays.Categorical
            'Ref' for Homozygous Reference
            'Het' for Heterozygous
            'Hom' for Homozygous Alt
            pd.NA for missing
            Raises an error if there is more than 1 alternate allele or ploidy is not 2
        """
        # TODO: Return multiple arrays for multiple alternate alleles?
        if len(self.variant.alleles) > 2:
            raise ValueError("Codominant encoding can only be used with one allele")
        if self.dtype.variant.ploidy != 2:
            raise ValueError(
                "Codominant encoding can only be used with diploid genotypes"
            )

        allele_sum = self.allele_idxs.sum(axis=1)
        categories = ["Ref", "Het", "Hom"]
        result = pd.Categorical(
            values=[categories[n] if n in {0, 1, 2} else None for n in allele_sum],
            categories=categories,
            ordered=True,
        )
        return result
