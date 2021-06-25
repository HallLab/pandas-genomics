"""
Test GenotypeArray Accessors
"""
import numpy as np
import pandas as pd
import pytest
from pandas._testing import (
    assert_series_equal,
)


def test_variant_score(data, data_for_encoding):
    assert pd.Series(data).genomics.variant.score == float(data.variant.score)
    assert pd.Series(data_for_encoding()).genomics.variant.score is None


def test_var_info(genotypearray_df):
    assert_series_equal(
        genotypearray_df.genomics.variant_info.iloc[0],
        genotypearray_df.iloc[:, 0].genomics.variant_info,
    )


def test_maf(data):
    assert pd.Series(data).genomics.maf == data.maf
    df = pd.DataFrame.from_dict(
        {n: pd.Series(data.copy()) for n in "ABC"}, orient="columns"
    )
    for colname in "ABC":
        df[colname].genomics.variant.id = colname
    df["D"] = np.ones(len(data))
    expected = pd.Series({"A": data.maf, "B": data.maf, "C": data.maf})
    assert_series_equal(df.genomics.maf, expected)


def test_hwe(data):
    assert pd.Series(data).genomics.hwe_pval == data.hwe_pval


@pytest.mark.parametrize(
    "filter_value, num_cols_left", [(None, 17), (0.05, 3), (0.10, 2)]
)
def test_filter_maf(genotypearray_df, filter_value, num_cols_left):
    if filter_value is None:
        result = genotypearray_df.genomics.filter_variants_maf()
    else:
        result = genotypearray_df.genomics.filter_variants_maf(filter_value)
    assert len(result.columns) == num_cols_left


@pytest.mark.parametrize(
    "filter_value, num_cols_left", [(None, 1), (0.05, 1), (1e-300, 2)]
)
def test_filter_hwe(ga_inhwe, ga_nothwe, filter_value, num_cols_left):
    data = pd.DataFrame({"yes": ga_inhwe, "no": ga_nothwe})
    data["num"] = [n for n in range(len(data))]
    if filter_value is None:
        result = data.genomics.filter_variants_hwe()
    else:
        result = data.genomics.filter_variants_hwe(filter_value)
    assert len(result.columns) == num_cols_left + 1
