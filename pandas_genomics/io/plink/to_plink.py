from typing import Optional

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

from pandas_genomics.arrays import GenotypeDtype


def to_plink(
    data: pd.DataFrame,
    output: str,
    phenotype_name: Optional[str] = None,
    phenotype_case: Optional[str] = None,
    phenotype_control: Optional[str] = None,
    id_prefix: str = "sample",
):
    """
    Save genetic data to plink v1 files (.bed, .bim, and .fam)

    Parameters
    ----------
    data: pd.DataFrame
        DataFrame containing GenotypeArrays to be saved.
    output: str
        Name to use for the output .bed, .bim, and .fam files
    phenotype_name: str, default None
        Optional column in `data` to be saved as the phenotype value in the .fam file.
    phenotype_case, phenotype_control
        String values indicating the category to be used as "case" or "control" for binary phenotypes.
        If provided, the phenotype must be categorical.
        If not provided, the phenotype is not encoded at all (assumed to be quantitative).
    id_prefix:
        If the data index is an integer index, this prefix will be added to generate IDs.

    Notes
    -----
    If the data index has the required columns (FID, IID, IID_father, IID_mother, sex, phenotype) the fam file will be created based on the index.
    If a phenotype name is provided, this will override any phenotype information in the index.
    If the data has a single index column this will be used (with the prefix) for FID and IID.  Defaults will be used for other .fam data

    """
    print(f"Saving genotype data to {output}")
    save_fam(
        data,
        output + ".fam",
        phenotype_name,
        phenotype_case,
        phenotype_control,
        id_prefix,
    )
    save_bim(data, output + ".bim")
    save_bed(data, output + ".bed")


def save_fam(
    data, output_fam, phenotype_name, phenotype_case, phenotype_control, id_prefix
):
    # Validate phenotype, if provided
    if phenotype_name is not None:
        if phenotype_name not in data.columns:
            raise ValueError(
                f"The phenotype ({phenotype_name}) was not found in the data"
            )
        else:
            phenotype_data = data[phenotype_name]

    # Check for full plink-style index, or create one
    if data.index.names == [
        "FID",
        "IID",
        "IID_father",
        "IID_mother",
        "sex",
        "phenotype",
    ]:
        # Recode sex
        fam_data = data.index.to_frame()
        fam_data["sex"].cat.rename_categories(
            {"male": 1, "female": 2, "unknown": 0}, inplace=True
        )
        # Update phenotype if provided
        if phenotype_name is not None:
            fam_data["phenotype"] = phenotype_data
    elif len(data.index.names) == 1:
        ids = pd.Series(data.index.values).apply(lambda i: f"{id_prefix}{i}")
        fam_data = pd.DataFrame.from_dict(
            {
                "FID": ids,
                "IID": ids,
                "IID_father": np.zeros(len(data)).astype(int),
                "IID_mother": np.zeros(len(data)).astype(int),
                "sex": np.zeros(len(data)).astype(int),
            }
        )
        if phenotype_name is None:
            fam_data["phenotype"] = np.ones(len(data)) * -9  # -9 is missing
        else:
            fam_data["phenotype"] = phenotype_data
    else:
        raise ValueError(
            "The data index must contain all 6 .fam file columns, or a single column"
        )

    # Optionally Encode phenotype
    if (phenotype_control is not None) and (phenotype_case is not None):
        if fam_data["phenotype"].dtype.name != "category":
            raise ValueError(
                "The phenotype must be categorical to utilize 'phenotype_control' and 'phenotype_case' parameters"
            )
        pheno_dict = {phenotype_control: 1, phenotype_case: 2}
        fam_data["phenotype"].cat.rename_categories(
            lambda c: pheno_dict.get(c, 0), inplace=True
        )

    fam_data.to_csv(output_fam, sep=" ", header=False, index=False)
    print(f"\tSaved information for {len(fam_data)} samples to {output_fam}")


def save_bim(data, output_bim):
    variants = [
        col_val.genomics.variant
        for col_name, col_val in data.iteritems()
        if GenotypeDtype.is_dtype(col_val.dtype)
    ]
    for var in variants:
        if len(var.alleles) != 2:
            raise ValueError(
                f"Variant {var.id} is not bialleleic (it has {len(var.alleles)} alleles) and therefore can't be saved in plink format."
            )
    var_dicts = [
        {
            "chromosome": var.chromosome,
            "variant_id": var.id,
            "position": 0,
            "coordinate": var.position,
            "allele1": var.alleles[1],  # alt
            "allele2": var.alleles[0],  # ref
        }
        for var in variants
    ]
    bim_data = pd.DataFrame(var_dicts)
    bim_data.to_csv(output_bim, sep="\t", header=False, index=False)
    print(f"\tSaved information for {len(bim_data)} variants to {output_bim}")


def save_bed(data, output_bed):
    # Get an array of bytes for each variant
    bytes = np.array(
        [
            gt_array_to_plink_bits(col_val)
            for col_name, col_val in data.iteritems()
            if GenotypeDtype.is_dtype(col_val.dtype)
        ]
    )
    # flatten into a single array
    bytes = bytes.flatten()
    # Add the first 3 bytes
    CORRECT_FIRST_BYTES = np.array([108, 27, 1], dtype="uint8")
    bytes = np.concatenate([CORRECT_FIRST_BYTES, bytes])
    # Write to file
    bytes.tofile(output_bed)
    print(f"\tSaved genotypes to {output_bed}")


def gt_array_to_plink_bits(gt_series):
    allele_ids = gt_series.array.allele_idxs
    # Replace missing with 0,1
    missing = gt_series.array.is_missing
    allele_ids[missing] = (0, 1)
    # Replace het with 1,0
    het_gt = allele_ids.sum(axis=1) == 1
    allele_ids[het_gt] = (1, 0)
    # Pad with zeros so it is divisible by 4
    pad_samples = len(allele_ids) % 4
    allele_ids = np.vstack([allele_ids, np.zeros(shape=(pad_samples, 2)).astype(int)])
    # Reshape into groups of 4 and flip each group
    allele_ids = np.flip(allele_ids.reshape((-1, 4, 2)), axis=1)
    # Flatten into a single array and pack bits into bytes
    allele_ids = np.packbits(allele_ids.flatten())
    return allele_ids


"""
    allele_ids = gt_array.genomics.allele_idxs
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
"""
