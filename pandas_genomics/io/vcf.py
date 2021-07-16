from pathlib import Path
from typing import Union

import pandas as pd
import numpy as np

from ..arrays import GenotypeArray, GenotypeDtype
from ..scalars import Variant, MISSING_IDX, Genotype


def from_vcf(
    filename: Union[str, Path], min_qual: float = 0, drop_filtered: bool = True
):
    """
    Load genetic data from a VCF or BCF file into a DataFrame

    Parameters
    ----------
    filename: str or Path
        vcf, vcf.gz, or bcf file.
    min_qual: float (default = 0)
        Skip loading variants with less than this quality
    drop_filtered: boolean (default = True)
        Skip loading variants with a FILTER value other than "PASS"

    Returns
    -------
    DataFrame
        Columns correspond to variants (named as {variant_number}_{variant ID}).
        Rows correspond to samples and index columns include sample information.

    Examples
    --------
    """
    from cyvcf2 import VCF  # Import here since installing htslib on Windows is tricky

    genotype_array_dict = dict()
    for var_num, vcf_variant in enumerate(VCF(filename)):  # or VCF('some.bcf')

        # Skip filtered variants unless drop_filtered is False
        if vcf_variant.FILTER is not None and drop_filtered:
            continue

        # Skip variants below the minimum quality
        if vcf_variant.QUAL < min_qual:
            continue

        if len(vcf_variant.ALT) >= MISSING_IDX:
            raise ValueError(
                f"Could not load {vcf_variant.ID} due to too many ALT alleles"
                f" ({len(vcf_variant.ALT)} > {MISSING_IDX-1})"
            )

        # Make variant
        variant = Variant(
            chromosome=vcf_variant.CHROM,
            position=vcf_variant.start,
            id=vcf_variant.ID,
            ref=vcf_variant.REF,
            alt=vcf_variant.ALT,
            ploidy=vcf_variant.ploidy,
            score=int(vcf_variant.QUAL),
        )
        dtype = GenotypeDtype(variant)

        # Collect genotypes
        allele_idxs = np.array(vcf_variant.genotypes)[:, :2]
        allele_idxs = np.where(allele_idxs == -1, MISSING_IDX, allele_idxs)
        gt_scores = vcf_variant.gt_quals
        # Convert genotype scores from float values to uint8 values
        gt_scores = np.where(gt_scores > 254, 254, gt_scores)  # Max Score
        gt_scores = np.where(gt_scores < 0, 255, gt_scores)  # Min Score (<0 is missing)
        gt_scores = np.where(gt_scores == -1, 255, gt_scores)  # Missing values
        gt_scores = gt_scores.round().astype("uint8")
        values = np.array(list(zip(allele_idxs, gt_scores)), dtype=dtype._record_type)

        # Make the GenotypeArray
        gt_array = GenotypeArray(values=values, dtype=dtype)
        # Make the variant name
        if gt_array.variant.id is None:
            var_name = f"Variant_{var_num}"
        else:
            var_name = gt_array.variant.id

        # Save to the dict
        genotype_array_dict[var_name] = gt_array

    df = pd.DataFrame.from_dict(genotype_array_dict)
    return df
