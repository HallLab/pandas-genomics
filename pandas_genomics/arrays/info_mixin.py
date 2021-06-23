from itertools import combinations_with_replacement

import numpy as np
from scipy.stats import chisquare

from pandas_genomics.arrays.utils import required_ploidy
from pandas_genomics.scalars import MISSING_IDX


class InfoMixin:
    """
    Genotype Mixin containing functions for calculating various information
    """

    @property
    def is_missing(self):
        """
        Boolean array: True if the sample is missing all alleles
        """
        return (self.allele_idxs == MISSING_IDX).all(axis=1)

    @property
    def is_homozygous(self):
        """
        Boolean array: True if the sample is homozygous for any allele
        """
        return (self.allele_idxs == self.allele_idxs[:, 0, np.newaxis]).all(axis=1)

    @property
    def is_heterozygous(self):
        """
        Boolean array: True if the sample is heterozygous for any alleles
        """
        return ~self.is_homozygous

    @property
    def is_homozygous_ref(self):
        """
        Boolean array: True if the sample is homozygous for the reference allele
        """
        return (self.allele_idxs == 0).all(axis=1)

    @property
    def is_homozygous_alt(self):
        """
        Boolean array: True if the sample is homozygous for any non-reference allele
        """
        return self.is_homozygous & ~self.is_homozygous_ref

    @property
    def maf(self) -> float:
        """
        Calculate the Minor Allele Frequency (MAF) for the most-frequent alternate allele.
        Missing alleles are ignored.
        """
        total_nonmissing_alleles = (self.allele_idxs != MISSING_IDX).sum().sum()
        if total_nonmissing_alleles == 0:
            # All genotypes missing
            return np.nan

        allele_counts = np.bincount(self.allele_idxs.flatten())
        if len(allele_counts) == 1:
            # All reference
            return 0.0
        elif len(allele_counts) == MISSING_IDX:
            # At least one missing allele is present: don't count them
            allele_counts = allele_counts[:-1]
        # Use highest alternate allele value
        return allele_counts[1:].max() / total_nonmissing_alleles

    @property
    @required_ploidy(2, np.nan)
    def hwe_pval(self) -> float:
        """
        Calculate the probability that the samples are in HWE for diploid variants

        Notes
        -----
        Generate expected counts using the allele frequencies and perform a chi-square test.
        Ignore any samples with missing alleles.
        Uses a typical number of degrees of freedom (the number of observed genotypes minus 1).
        Returns np.nan if any expected counts are < 5
        """
        # Take nonmissing allele indexes
        nonmissing_aidxs = self.allele_idxs[self.allele_idxs.max(axis=1) != MISSING_IDX]
        if len(nonmissing_aidxs) == 0:
            return np.nan

        # Get allele counts and frequency
        allele_counts = np.bincount(nonmissing_aidxs.flatten())
        total_gt = len(nonmissing_aidxs)
        total_alleles = total_gt * 2
        if total_gt < 2:
            return np.nan  # Too few samples to calculate
        if len(allele_counts) == 1:
            return 1.0  # All Reference
        allele_freqs = allele_counts / (total_alleles)

        # Chisq test
        expected = []
        observed = []
        for a1, a2 in combinations_with_replacement(range(len(allele_freqs)), 2):
            a1_freq = allele_freqs[a1]
            a2_freq = allele_freqs[a2]
            if a1 != a2:
                # Heterozygous
                expected.append(int(a1_freq * a2_freq * total_gt * 2))
            else:
                # Homozygous
                expected.append(int(a1_freq * a2_freq * total_gt))
            observed.append((nonmissing_aidxs == (a1, a2)).all(axis=1).sum())
        # Return NaN if any expected counts are < 5
        if min(expected) < 5:
            return np.nan
        chi, pval = chisquare(observed, expected)
        return pval
