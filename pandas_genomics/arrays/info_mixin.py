from itertools import combinations_with_replacement

import numpy as np
import pandas as pd
from scipy.stats import chisquare

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

    @property
    def hwe_pval(self) -> float:
        """
        Calculate the probability that the samples are in HWE.

        Generate expected counts using the allele frequencies and perform a chi-square test.
        Ignore any samples with missing alleles.
        Uses a typical number of degrees of freedom (the number of observed genotypes minus 1)
        """
        # Only calculate for diploid, otherwise the math is a lot more complicated.
        # Potentially support this in the future
        if self.variant.ploidy != 2:
            return np.nan

        # Take nonmissing allele indexes
        nonmissing_aidxs = self.allele_idxs[self.allele_idxs.max(axis=1) != MISSING_IDX]
        if len(nonmissing_aidxs) == 0:
            return np.nan

        # Get allele counts and frequency
        allele_counts = np.bincount(nonmissing_aidxs.flatten())
        total_gt = len(nonmissing_aidxs)
        total_alleles = total_gt * 2
        if len(allele_counts) == 1:
            return 1.0  # All Reference
        if total_gt < 2:
            return np.nan  # Too few samples to calculate
        allele_freqs = allele_counts / (total_alleles)

        # Chisq test
        expected = []
        observed = []
        for a1, a2 in combinations_with_replacement(range(len(allele_freqs)), 2):
            a1_freq = allele_freqs[a1]
            a2_freq = allele_freqs[a2]
            expected.append(a1_freq * a2_freq * total_alleles)
            observed.append(
                (
                    (nonmissing_aidxs == (a1, a2)).all(axis=1)
                    | (nonmissing_aidxs == (a1, a2)).all(axis=1)
                ).sum()
            )
        chi, pval = chisquare(observed, expected)
        return pval
