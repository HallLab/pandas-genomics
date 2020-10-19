from pathlib import Path
from typing import Union

import pandas as pd

from ..arrays import GenotypeArray
from ..scalars import Variant


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
        # TODO: Should FILTER or QUAL be stored in the GenotypeArray?

        # Skip filtered variants unless drop_filtered is True
        if vcf_variant.FILTER is not None and drop_filtered:
            continue

        # Skip variants below the minimum quality
        if vcf_variant.QUAL < min_qual:
            continue

        # Make variant
        variant = Variant(
            chromosome=vcf_variant.CHROM,
            position=vcf_variant.start,
            id=vcf_variant.ID,
            ref=vcf_variant.REF,
            alt=vcf_variant.ALT,
        )
        # Make the GenotypeArray
        gt_array = GenotypeArray(
            values=[
                variant.make_genotype_from_vcf_record(vcf_record)
                for vcf_record in vcf_variant.genotypes
            ]
        )
        # Make the variant name
        if gt_array.variant.id is None:
            var_name = f"Variant_{var_num}"
        else:
            var_name = gt_array.variant.id

        # Save to the dict
        genotype_array_dict[var_name] = gt_array

    df = pd.DataFrame.from_dict(genotype_array_dict)
    return df
