from typing import Optional

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

from pandas_genomics.arrays import GenotypeDtype


def to_plink(
        data: pd.DataFrame,
        output: str,
        phenotype: Optional[str] = None,
):
    """
    Save genetic data to plink v1 files (.bed, .bim, and .fam)

    Parameters
    ----------
    data: pd.DataFrame
        DataFrame containing GenotypeArrays to be saved.
    output: str
        Name to use for the output .bed, .bim, and .fam files
    phenotype: str, default None
        Optional column in `data` to be saved as the phenotype value in the .fam file.

    Notes
    -----
    If the data index has the required columns (FID, IID, IID_father, IID_mother, sex, phenotype) the fam file will
      be created based on the index.
    If a phenotype name is provided, this will override any phenotype information in the index.
    If the data has a single index column this will be used for IID and defaults will be used for other .fam data

    """
    save_fam(data, phenotype, output+".fam")
    save_bim(data, output+".bim")
    save_bed(data, output+".bed")


def save_fam(data, phenotype, output_fam):
    # Validate phenotype, if provided
    if phenotype is not None:
        if phenotype not in data.columns:
            raise ValueError(f"The phenotype ({phenotype}) was not found in the data")
        elif not is_numeric_dtype(data[phenotype]):
            raise ValueError("Only numeric phenotypes are allowed in .fam files.")
        else:
            phenotype = data["phenotype"]

    # Check for full plink-style index, or create one
    if data.index.names == ["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"]:
        # Recode sex
        fam_data = data.index.to_frame()
        fam_data["sex"].cat.rename_categories({"male": 1, "female": 2, "unknown": 0}, inplace=True)
        # Update phenotype if provided
        if phenotype is not None:
            fam_data["phenotype"] = phenotype
    elif len(data.index.names) == 1:
        if 0 in data.index:
            raise ValueError("Can't use the index for .fam file IID because it contains '0'")
        fam_data = pd.DataFrame.from_dict({
            "FID": data.index,
            "IID": data.index,
            "IID_father": np.zeros(len(data)),
            "IID_mother": np.zeros(len(data)),
            "sex": np.zeros(len(data))})
        if phenotype is None:
            fam_data["phenotype"] = np.ones(len(data)) * -9  # -9 is missing
        else:
            fam_data["phenotype"] = phenotype
    else:
        raise ValueError("The data index must contain all 6 .fam file columns, or a single column")

    fam_data.to_csv(output_fam, sep=" ", header=False, index=False)


def save_bim(data, output_bim):
    variants = [col_val.genomics.variant for col_name, col_val in data.iteritems() if GenotypeDtype.is_dtype(col_val.dtype)]
    for var in variants:
        if len(var.alleles) != 2:
            raise ValueError(f"Variant {var.id} is not bialleleic (it has {len(var.alleles)} alleles) and therefore can't be saved in plink format.")
    var_dicts = [{
        "chromosome": var.chromosome,
        "variant_id": var.id,
        "position": 0,
        "coordinate": var.position,
        "allele1": var.alleles[1],  # alt
        "allele2": var.alleles[0]  # ref
        } for var in variants]
    bim_data = pd.DataFrame(var_dicts)
    bim_data.to_csv(output_bim, sep="\t", header=False, index=False)

def save_bed(data, output_bed):
    pass
