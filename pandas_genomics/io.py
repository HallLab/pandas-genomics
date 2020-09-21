from pathlib import Path
import pandas as pd
from .dtypes.genotype import GenotypeDtype, GenotypeArray, Variant


def from_plink(bed_file: str):
    """
    Load genetic data from plink files (.bed, .bim, and .fam) into a DataFrame

    Parameters
    ----------
    bed_file: str or Path
        PLINK .bed file.  .bim and .fam files with the same name and location must also exist.

    Returns
    -------
    DataFrame
        Index columns include sample information.
        Columns correspond to variants and rows correspond to samples.

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
    df['sex'] = df['sex'].astype('category')
    df['sex'] = df['sex'].cat.rename_categories({1: 'male', 2: 'female', 0: 'unknown'})
    # 'Phenotype' is not encoded because it may be more complicated than just case/control
    num_samples = len(df)
    print(f"\tLoaded information for {num_samples} samples from '{fam_file.name}'")

    # Load bim file (PLINK extended MAP file)
    variant_info = pd.read_table(bim_file, header=None, sep="\t")
    variant_info.columns = ["chromosome", "variant_id", "position", "coordinate", "allele1", "allele2"]
    # chromosome is a category
    variant_info['chromosome'] = variant_info['chromosome'].astype('category')
    num_variants = len(variant_info)
    print(f"\tLoaded information for {num_variants} variants from '{bim_file.name}'")

    # Load bed file (PLINK binary biallelic genotype table) and add info to the df
    CORRECT_FIRST_BYTES = b"\x6c\x1b\x01"
    with bed_file.open('rb') as f:
        first3bytes = f.read(3)  # Read first three bytes and confirm they are correct
        if first3bytes != CORRECT_FIRST_BYTES:
            raise ValueError(f"The first 3 bytes {bed_file.name} were not correct.  The file may be corrupted.")
        # Determine chunk size
        chunk_size = num_samples // 4
        if num_samples % 4 > 0:
            chunk_size += 1
        # Read through the file
        for v_idx in range(num_variants):
            variant_info_dict = variant_info.iloc[v_idx].to_dict()
            variant_id = str(variant_info_dict['variant_id'])
            variant = Variant(variant_id=variant_id,
                              chromosome=str(variant_info_dict['chromosome']),
                              coordinate=int(variant_info_dict['coordinate']),
                              alleles=[str(variant_info_dict['allele1']), str(variant_info_dict['allele2'])])
            genotypes = []
            chunk = f.read(chunk_size)  # Encoded chunk of results for each variant
            for byte in chunk:
                # for each byte, get 2 bits at a time in reverse order (as a string, so '00', '01', '10', or '11')
                bitstrings = [f"{byte:08b}"[i:i+2] for i in range(0, 8, 2)][::-1]
                genotypes.extend([variant.make_genotype_from_plink_bits(bs) for bs in bitstrings])
            # Remove nonexistent samples at the end
            genotypes = genotypes[:num_samples]
            df[variant_id] = GenotypeArray(values=genotypes, dtype=GenotypeDtype(variant))
    print(f"\tLoaded genotypes from '{bed_file.name}'")

    # Set sample info as the index
    df = df.set_index(["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"])

    return df
