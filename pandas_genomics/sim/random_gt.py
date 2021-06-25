from typing import List, Union

import numpy as np

from pandas_genomics.arrays import GenotypeArray, GenotypeDtype
from pandas_genomics.scalars import Variant


def generate_random_gt(
    variant: Variant,
    alt_allele_freq: Union[List[float], float],
    n: int = 1000,
    random_seed: int = 1855,
) -> GenotypeArray:
    """
    Simulate random genotypes according to the provided allele frequencies

    Parameters
    ----------
    variant: Variant
    alt_allele_freq: float or List[float]
      Allele frequencies for each alternate allele in the variant (Bialleleic variants may specify a single float value)
    n: int, default 1000
        How many genotypes to simulate
    random_seed: int, default 1855

    Returns
    -------
    GenotypeArray

    """
    # Validate frequencies
    if isinstance(alt_allele_freq, float):
        # Convert it into a list
        alt_allele_freq = [
            alt_allele_freq,
        ]
    if len(alt_allele_freq) != len(variant.alleles) - 1:
        raise ValueError(
            f"The number of provided frequencies ({len(alt_allele_freq)}) doesn't match"
            f" the number of alternate alleles in the variant ({len(variant.alleles)-1})."
        )
    if sum(alt_allele_freq) > 1.0:
        raise ValueError(
            f"The provided frequencies must not sum to > 1.0 (sum was {sum(alt_allele_freq):.3e})"
        )

    # Set remaining odds to the reference allele
    allele_freq = [
        1 - sum(alt_allele_freq),
    ] + alt_allele_freq

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
