from pathlib import Path

from cyvcf2 import VCF
import pandas as pd

from ..arrays import GenotypeDtype, GenotypeArray
from ..scalars import Variant


def from_vcf(filename: str, min_qual: float = 0, drop_filtered: bool = True):
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
    genotype_array_list = []
    for vcf_variant in VCF(filename):  # or VCF('some.bcf')
        # TODO: Should FILTER or QUAL be stored in the GenotypeArray?

        # Skip filtered variants unless drop_filtered is True
        if vcf_variant.FILTER is not None and drop_filtered:
            continue

        # Skip variants below the minimum quality
        if vcf_variant.QUAL < min_qual:
            continue

        variant = Variant(chromosome = vcf_variant.CHROM,
                          position = vcf_variant.start,
                          id=vcf_variant.ID,
                          ref=vcf_variant.REF,
                          alt=vcf_variant.ALT)

        alleles = vcf_variant.gt_bases
        print()
        