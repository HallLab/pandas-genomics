from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np
from ..arrays import GenotypeDtype, GenotypeArray
from ..scalars import Variant, MISSING_IDX, Genotype


def from_plink(
    bed_file: str, swap_alleles: bool = False, max_variants: Optional[int] = None
):
    """
    Load genetic data from plink files (.bed, .bim, and .fam) into a DataFrame.

    Parameters
    ----------
    bed_file: str or Path
        PLINK .bed file.  .bim and .fam files with the same name and location must also exist.
    swap_alleles: bool
        False by default, in which case "allele2" (usually major) in the bim file is considered the "reference" allele.
        If True, "allele1" (usually minor) is considered the "reference" allele.
    max_variants: Optional[int]
        If provided, only load this number of variants

    Returns
    -------
    DataFrame
        Columns correspond to variants (named as {variant_number}_{variant ID}).
        Rows correspond to samples and index columns include sample information.

    Notes
    -----
    Plink files encode all variants as diploid (2n) and utilize "missing" alleles if the variant is actually haploid

    Examples
    --------
    """
    # All three files are used
    bed_file = Path(bed_file)
    name = bed_file.stem
    bim_file = bed_file.parent / f"{name}.bim"
    fam_file = bed_file.parent / f"{name}.fam"

    # Make sure each file exists
    if not Path(bed_file).exists():
        raise ValueError(f"The .bed file was not found\n\t{str(bed_file)}")
    if not Path(bim_file).exists():
        raise ValueError(f"The .bim file was not found\n\t{str(bim_file)}")
    if not Path(fam_file).exists():
        raise ValueError(f"The .fam file was not found\n\t{str(fam_file)}")

    print(f"Loading genetic data from '{bed_file.stem}'")

    # Load fam file (PLINK sample information file)
    df = pd.read_table(fam_file, header=None, sep=" ")
    df.columns = ["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"]
    # Update 'sex'
    df["sex"] = df["sex"].astype("category")
    df["sex"] = df["sex"].cat.rename_categories({1: "male", 2: "female", 0: "unknown"})
    # 'Phenotype' is not encoded because it may be more complicated than just case/control
    num_samples = len(df)
    print(f"\tLoaded information for {num_samples} samples from '{fam_file.name}'")

    # Load bim file (PLINK extended MAP file)
    variant_info = pd.read_table(bim_file, header=None, sep="\t")
    # Note 'position' is in centimorgans, 'coordinate' is what pandas-genomics refers to as 'position' (in base-pairs)
    variant_info.columns = [
        "chromosome",
        "variant_id",
        "position",
        "coordinate",
        "allele1",
        "allele2",
    ]
    # chromosome is a category
    variant_info["chromosome"] = variant_info["chromosome"].astype("category")
    num_variants = len(variant_info)

    # Limit num_variants
    if max_variants is not None:
        if max_variants < 1:
            raise ValueError(f"'max_variants' set to an invalid value: {max_variants}")
        else:
            num_variants = max_variants

    print(f"\tLoaded information for {num_variants} variants from '{bim_file.name}'")

    # Load bed file (PLINK binary biallelic genotype table)
    gt_bytes = np.fromfile(bed_file, dtype="uint8")
    # Ensure the file is valid
    CORRECT_FIRST_BYTES = np.array([108, 27, 1], dtype="uint8")
    if not (gt_bytes[:3] == CORRECT_FIRST_BYTES).all():
        raise ValueError(
            f"The first 3 bytes {bed_file.name} were not correct.  The file may be corrupted."
        )
    gt_bytes = gt_bytes[3:]
    # Divide array into one row per variant
    # (have to reshape using num_samples since num_variants may be set lower)
    chunk_size = num_samples // 4
    if num_samples % 4 > 0:
        chunk_size += 1
    gt_bytes = gt_bytes.reshape(-1, chunk_size)
    # Process each variant
    for v_idx in range(num_variants):
        variant_info_dict = variant_info.iloc[v_idx].to_dict()
        variant_id = str(variant_info_dict["variant_id"])
        a1 = str(variant_info_dict["allele1"])
        a2 = str(variant_info_dict["allele2"])
        # 0 indicates a missing allele
        if a2 == "0":
            ref = None
        if a1 == "0":
            a1 = None
        else:
            a1 = [a1]  # pass as list
        variant = Variant(
            chromosome=str(variant_info_dict["chromosome"]),
            position=int(variant_info_dict["coordinate"]),
            id=variant_id,
            ref=a2,
            alt=a1,
            ploidy=2,
        )
        # Each byte (8 bits) is a concatenation of two bits per sample for 4 samples
        # These are ordered from right to left, like (sample4, sample3, sample2, sample1)
        # Convert each byte into 4 2-bits and flip them to order samples correctly
        genotypes = np.flip(np.unpackbits(gt_bytes[v_idx]).reshape(-1, 4, 2), axis=1)
        # flatten the middle dimension to give a big list of genotypes in the correct order and
        # remove excess genotypes at the end that are padding rather than real samples
        genotypes = genotypes.reshape(-1, 2)[:num_samples]
        # Replace 0, 1 with missing (1, 0 is heterozygous)
        missing_gt = (genotypes == (0, 1)).all(axis=1)
        genotypes[missing_gt] = (MISSING_IDX, MISSING_IDX)
        # Replace 1, 0 with 0, 1 for heterozygous so the reference allele is first
        het_gt = (genotypes == (1, 0)).all(axis=1)
        genotypes[het_gt] = (0, 1)
        # Create GenotypeArray representation of the data
        dtype = GenotypeDtype(variant)
        scores = np.empty(num_samples)
        scores[:] = np.nan
        data = np.array(list(zip(genotypes, scores)), dtype=dtype._record_type)
        gt_array = GenotypeArray(values=data, dtype=dtype)
        if swap_alleles:
            gt_array.set_reference(1)
        df[f"{v_idx}_{variant_id}"] = gt_array
    print(f"\tLoaded genotypes from '{bed_file.name}'")

    # Set sample info as the index
    df = df.set_index(["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"])

    return df
