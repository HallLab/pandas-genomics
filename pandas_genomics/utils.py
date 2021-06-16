from typing import Union, Optional, List
import pandas as pd
import statsmodels.api as sm

from .arrays import GenotypeDtype


def generate_weighted_encodings(
    genotypes: Union[pd.Series, pd.DataFrame],
    data: pd.DataFrame,
    outcome_variable: str,
    covariates: Optional[List[str]] = None,
):
    """
    Calculate alpha values to be used in weighted encoding

    Parameters
    ----------
    genotypes:
        A GenotypeArray Series or DataFrame
    data:
        Data to be used in the regression, including the outcome and covariates
    outcome_variable:
        The variable to be used as the output (y) of the regression
    covariates:
        Other variables to be included in the regression formula

    Returns
    -------
    Dict
      Variant ID: str
      Alpha Value - used for heterozygous genotypes
      Ref Allele - which allele is considered reference
      Alt Allele - which allele is considered alternate
      Minor Allele Frequency - MAF of data used during calculation of alpha values

    Notes
    -----
    See [1]_ for more information about weighted encoding.

    References
    ----------
    .. [1] Hall, Molly A., et al.
           "Novel EDGE encoding method enhances ability to identify genetic interactions."
           PLoS genetics 17.6 (2021): e1009534.
    """
    # Validate parameters
    if covariates is None:
        covariates = []
    # Covariates must be a list
    if type(covariates) != list:
        raise ValueError("'covariates' must be specified as a list or set to None")

    # Extract specific data
    try:
        data = data[
            [
                outcome_variable,
            ]
            + covariates
        ]
    except KeyError as e:
        raise ValueError(f"Missing variable in provided data: {e}")

    # Check Types to determine which kind of regression to run
    dtypes = _get_types(data)

    outcome_type = dtypes.get(outcome_variable)
    if outcome_type == "continuous":
        family = sm.families.Gaussian(link=sm.families.links.identity())
        use_t = True
    elif outcome_type == "binary":
        # Use the order according to the categorical
        counts = data[outcome_variable].value_counts().to_dict()
        categories = data[outcome_variable].cat.categories
        codes, categories = zip(*enumerate(categories))
        data[outcome_variable].replace(categories, codes, inplace=True)
        print(
            f"Binary Outcome (family = Binomial): '{outcome_variable}'\n"
            f"\t{counts[categories[0]]:,} occurrences of '{categories[0]}' coded as 0\n"
            f"\t{counts[categories[1]]:,} occurrences of '{categories[1]}' coded as 1"
        )
        family = sm.families.Binomial(link=sm.families.links.logit())
        use_t = False

    # Check for missing outcomes
    na_outcome_count = data[outcome_variable].isna().sum()
    if na_outcome_count > 0:
        raise ValueError(f"{na_outcome_count} samples are missing an outcome value")

    # Ensure genotypes data is actually all genotypes
    if not genotypes.dtypes.apply(lambda dt: GenotypeDtype.is_dtype(dt)).all():
        incorrect = genotypes.dtypes[
            ~genotypes.dtypes.apply(lambda dt: GenotypeDtype.is_dtype(dt))
        ]
        raise AttributeError(
            f"Incompatible datatypes: all columns must be a GenotypeDtype: {incorrect}"
        )

    # Merge genotypes and data
    for col in list(data):
        if col in list(genotypes):
            raise ValueError(
                "Outcome and covariate names should not exist in `genotypes`: Check '{col}'"
            )
    gt_col_names = list(genotypes)
    merged = genotypes.merge(data, how="inner", left_index=True, right_index=True)
    if len(merged) == 0:
        raise ValueError(
            "Unable to merge the genotypes with the data.  Check the index values."
        )
    elif len(merged) < len(genotypes):
        raise ValueError(
            f"Only {len(merged):,} of {len(genotypes):,} genotypes were merged to the data.  Check the index values."
        )

    # Run regressions
    for gt in gt_col_names:
        encoded = merged[gt].genomics.encode_codominant()
        df = pd.concat([data, encoded], axis=1)
        formula = f"Q('{outcome_variable}') ~ 1"
        if len(covariates) > 0:
            formula += " + "
            formula += " + ".join([f"Q('{v}')" for v in varying_covars])


def _get_types(data: pd.DataFrame):
    """Categorize data into 'unknown', 'constant', 'binary', 'categorical', or 'continuous'"""
    dtypes = pd.Series("unknown", index=data.columns)

    # Set binary and categorical
    data_catbin = data.loc[:, data.dtypes == "category"]
    if len(data_catbin.columns) > 0:
        # Constant
        constant_cols = data_catbin.apply(lambda col: len(col.cat.categories) == 1)
        constant_cols = constant_cols[constant_cols].index
        dtypes.loc[constant_cols] = "constant"
        # Binary
        bin_cols = data_catbin.apply(lambda col: len(col.cat.categories) == 2)
        bin_cols = bin_cols[bin_cols].index
        dtypes.loc[bin_cols] = "binary"
        # Categorical
        cat_cols = data_catbin.apply(lambda col: len(col.cat.categories) > 2)
        cat_cols = cat_cols[cat_cols].index
        dtypes.loc[cat_cols] = "categorical"

    # Set continuous
    cont_cols = data.dtypes.apply(lambda dt: pd.api.types.is_numeric_dtype(dt))
    cont_cols = cont_cols[cont_cols].index
    dtypes.loc[cont_cols] = "continuous"

    return dtypes
