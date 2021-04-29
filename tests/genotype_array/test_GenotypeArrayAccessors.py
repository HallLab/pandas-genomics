"""
Test GenotypeArray Accessors
"""

import pandas as pd
import pytest
from pandas._testing import (
    assert_series_equal,
)


def test_variant_score(data, data_for_encoding):
    assert pd.Series(data).genomics.variant.score == float(data.variant.score)
    assert pd.Series(data_for_encoding).genomics.variant.score is None


def test_var_info(genotypearray_df):
    assert_series_equal(
        genotypearray_df.genomics.variant_info.iloc[0],
        genotypearray_df.iloc[:, 0].genomics.variant_info,
    )


def test_maf(data):
    assert pd.Series(data).genomics.maf == data.maf
    df = pd.DataFrame.from_dict({n: pd.Series(data) for n in "ABC"}, orient="columns")
    expected = pd.Series({n: data.maf for n in "ABC"})
    assert_series_equal(df.genomics.maf, expected)


def test_hwe(data):
    assert pd.Series(data).genomics.hwe_pval == data.hwe_pval


@pytest.mark.parametrize(
    "filter_value, num_vars_left", [(None, 15), (0.05, 1), (0.10, 0)]
)
def test_filter_maf(genotypearray_df, filter_value, num_vars_left):
    if filter_value is None:
        result = genotypearray_df.genomics.filter_variants_maf()
    else:
        result = genotypearray_df.genomics.filter_variants_maf(filter_value)
    assert len(result.columns) == num_vars_left


@pytest.mark.parametrize(
    "filter_value, num_vars_left", [(None, 1), (0.04, 1), (1e-18, 17)]
)
def test_filter_maf(genotypearray_df, filter_value, num_vars_left):
    if filter_value is None:
        result = genotypearray_df.genomics.filter_variants_hwe()
    else:
        result = genotypearray_df.genomics.filter_variants_hwe(filter_value)
    assert len(result.columns) == num_vars_left
