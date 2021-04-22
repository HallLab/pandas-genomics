from typing import List

import numpy as np

from pandas_genomics.arrays import GenotypeArray, GenotypeDtype
from pandas_genomics.scalars import Variant


def generate_random_gt(
    variant: Variant, allele_freq: List[float], n: int = 1000, random_seed: int = 1855
) -> GenotypeArray:
    """
    Simulate random genotypes according to the provided allele frequencies

    Parameters
    ----------
    variant: Variant
    allele_freq: List[float]
      Allele frequencies for each allele in the variant
    n: int, default 1000
        How many genotypes to simulate
    random_seed: int, default 1855

    Returns
    -------
    GenotypeArray

    """
    # Validate frequencies
    if len(allele_freq) != len(variant.alleles):
        raise ValueError(
            f"The number of provided frequencies ({len(allele_freq)}) doesn't match"
            f" the number of alleles in the variant ({len(variant.alleles)})."
        )
    if sum(allele_freq) != 1.0:
        raise ValueError(
            f"The provided frequencies must add up to 1.0 (sum was {sum(allele_freq):.3f})"
        )

    # Choose gts
    np.random.seed(random_seed)
    genotypes = np.random.choice(
        range(len(variant.alleles)), p=allele_freq, size=(n, variant.ploidy)
    )

    # Create GenotypeArray representation of the data
    dtype = GenotypeDtype(variant)
    scores = np.empty(n)
    scores[:] = np.nan
    data = np.array(list(zip(genotypes, scores)), dtype=dtype._record_type)
    gt_array = GenotypeArray(values=data, dtype=dtype)

    return gt_array
