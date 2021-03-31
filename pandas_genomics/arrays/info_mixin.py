import numpy as np
import pandas as pd

from pandas_genomics.scalars import MISSING_IDX


class InfoMixin:
    """
    Genotype Mixin containing functions for calculating various information
    """

    @property
    def is_missing(self):
        """
        Get a boolean array indicating which samples are missing
        """
        return np.array((self.allele_idxs == (MISSING_IDX, MISSING_IDX)).all(axis=1))

    @property
    def maf(self) -> float:
        """
        Calculate the Minor Allele Frequency (MAF) for the most-frequent alternate allele.
        Missing genotypes are ignored.
        """
        total_alleles = (1 - self.is_missing).sum() * self.variant.ploidy
        if total_alleles == 0:
            # All genotypes missing
            return np.nan

        allele_counts = np.bincount(self.allele_idxs.flatten())
        if len(allele_counts) == 1:
            # All reference
            return 0.0
        elif len(allele_counts) == MISSING_IDX:
            # At least one missing allele is present: don't count it
            allele_counts = allele_counts[:-1]
        # Use highest alternate allele value
        return allele_counts[1:].max() / total_alleles
