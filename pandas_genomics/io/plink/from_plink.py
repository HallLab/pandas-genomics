from pathlib import Path
from typing import Optional, Union

import pandas as pd
import numpy as np
from ...arrays import GenotypeDtype, GenotypeArray
from ...scalars import Variant, MISSING_IDX


def from_plink(
    input: Union[str, Path],
    swap_alleles: bool = False,
    max_variants: Optional[int] = None,
    categorical_phenotype: bool = True,
):
    """
    Load genetic data from plink v1 files (.bed, .bim, and .fam) into a DataFrame.

    Parameters
    ----------
    input: str or Path
        PLINK sample (no extension):  .bed, .bim and .fam files with the same name and location must exist.
    swap_alleles: bool
        False by default, in which case "allele2" (usually major) in the bim file is considered the "reference" allele.
        If True, "allele1" (usually minor) is considered the "reference" allele.
    max_variants: Optional[int]
        If provided, only load this number of variants
    categorical_phenotype: bool, True by default
        If True, the phenotype is encoded as a categorical when loaded (1 = "Control", 2 = "Case", otherwise missing.
        If False, load values directly.

    Returns
    -------
    DataFrame
        Columns correspond to variants (named as {variant_number}_{variant ID}).
        Rows correspond to samples and index columns include sample information.

    Notes
    -----
    Plink v1 files encode all variants as diploid (2n) and utilize "missing" alleles if the variant is actually haploid

    Examples
    --------
    """
    # All three files are used
    input = str(input)  # coerce to string in order to add extensions
    bed_file = Path(input + ".bed")
    bim_file = Path(input + ".bim")
    fam_file = Path(input + ".fam")

    # Make sure each file exists
    if not Path(bed_file).exists():
        raise ValueError(f"The .bed file was not found\n\t{str(bed_file)}")
    if not Path(bim_file).exists():
        raise ValueError(f"The .bim file was not found\n\t{str(bim_file)}")
    if not Path(fam_file).exists():
        raise ValueError(f"The .fam file was not found\n\t{str(fam_file)}")

    print(f"Loading genetic data from '{bed_file.stem}'")

    # Load fam file
    df = load_sample_info(fam_file, categorical_phenotype)
    # Lod bim file
    variant_list = load_variant_info(bim_file, max_variants)  # Load bed file
    gt_array_dict = load_genotypes(
        bed_file, variant_list, num_samples=len(df), swap_alleles=swap_alleles
    )

    # Merge with sample allele index
    df = pd.concat([df, pd.DataFrame.from_dict(gt_array_dict)], axis=1)
    df = df.set_index(["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"])

    return df


def load_sample_info(fam_file, categorical_phenotype):
    """Load fam file (PLINK sample information file) into a df"""
    df = pd.read_table(fam_file, header=None, sep=" ")
    df.columns = ["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"]
    # Update 'sex'
    df["sex"] = df["sex"].astype("category")
    df["sex"] = df["sex"].cat.rename_categories({1: "male", 2: "female", 0: "unknown"})
    # Encode the phenotype
    DEFAULT_CAT_MAP = {1: "Control", 2: "Case"}
    if categorical_phenotype:
        df["phenotype"] = df["phenotype"].astype("category")
        df["phenotype"].cat.rename_categories(DEFAULT_CAT_MAP, inplace=True)
        df.loc[~df["phenotype"].isin(DEFAULT_CAT_MAP.values()), "phenotype"] = None
    print(f"\tLoaded information for {len(df)} samples from '{fam_file.name}'")
    return df


def load_variant_info(bim_file, max_variants):
    """Load bim file (PLINK extended MAP file) into a list of variants"""
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
    # Limit num_variants
    if max_variants is not None:
        if max_variants < 1:
            raise ValueError(f"'max_variants' set to an invalid value: {max_variants}")
        else:
            variant_info = variant_info.iloc[:max_variants]
    variant_list = [create_variant(row) for idx, row in variant_info.iterrows()]
    print(
        f"\tLoaded information for {len(variant_list)} variants from '{bim_file.name}'"
    )
    return variant_list


def create_variant(variant_info_row):
    variant_id = str(variant_info_row["variant_id"])
    a1 = str(variant_info_row["allele1"])
    a2 = str(variant_info_row["allele2"])
    # 0 indicates a missing allele
    if a2 == "0":
        a2 = None
    if a1 == "0":
        a1 = None
    else:
        a1 = [a1]  # pass as list
    # Ensure chromosome is None instead of nan
    if np.isnan(variant_info_row["chromosome"]):
        chromosome = None
    else:
        chromosome = str(variant_info_row["chromosome"])
    variant = Variant(
        chromosome=chromosome,
        position=int(variant_info_row["coordinate"]),
        id=variant_id,
        ref=a2,
        alt=a1,
        ploidy=2,
    )
    return variant


def load_genotypes(bed_file, variant_list, num_samples, swap_alleles):
    """Load bed file (PLINK binary biallelic genotype table) into a dictionary of name:GenotypeArray"""
    gt_bytes = np.fromfile(bed_file, dtype="uint8")
    # Ensure the file is valid
    CORRECT_FIRST_BYTES = np.array([108, 27, 1], dtype="uint8")
    if not (gt_bytes[:3] == CORRECT_FIRST_BYTES).all():
        raise ValueError(
            f"The first 3 bytes {bed_file.name} were not correct.  The file may be corrupted."
        )
    gt_bytes = gt_bytes[3:]
    # Divide array into one row per variant
    chunk_size = num_samples // 4
    if num_samples % 4 > 0:
        chunk_size += 1
    gt_bytes = gt_bytes.reshape(-1, chunk_size)
    # Process each variant
    gt_array_dict = {}
    for v_idx, variant in enumerate(variant_list):
        variant_gt_bytes = gt_bytes[v_idx]
        gt_array = create_gt_array(num_samples, variant_gt_bytes, variant)
        if swap_alleles:
            gt_array.set_reference(1)
        gt_array_dict[f"{v_idx}_{gt_array.variant.id}"] = gt_array
    print(f"\tLoaded genotypes from '{bed_file.name}'")
    return gt_array_dict


def create_gt_array(num_samples, variant_gt_bytes, variant):
    # Each byte (8 bits) is a concatenation of two bits per sample for 4 samples
    # These are ordered from right to left, like (sample4, sample3, sample2, sample1)
    # Convert each byte into 4 2-bits and flip them to order samples correctly
    genotypes = np.flip(np.unpackbits(variant_gt_bytes).reshape(-1, 4, 2), axis=1)
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
    scores = np.ones(num_samples) * MISSING_IDX  # Missing Scores
    data = np.array(list(zip(genotypes, scores)), dtype=dtype._record_type)
    gt_array = GenotypeArray(values=data, dtype=dtype)
    return gt_array
