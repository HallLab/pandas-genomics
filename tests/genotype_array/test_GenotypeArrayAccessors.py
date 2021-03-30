"""
Test GenotypeArray methods
"""

import pandas as pd
import pytest
from pandas._testing import assert_series_equal, assert_extension_array_equal


def test_variant_score(data, data_for_encoding):
    assert pd.Series(data).genotype.variant.score == 30.0
    assert pd.Series(data_for_encoding).genotype.variant.score is None


@pytest.mark.xfail(raises=ValueError)
def test_encoding_extra_alt(data):
    """Fail encoding when there are multiple alt alleles"""
    data.encode_additive()


def test_encoding_additive(data_for_encoding):
    # Test arrays directly
    expected = pd.array([0, 1, 2, None], dtype="UInt8")
    result = data_for_encoding.encode_additive()
    assert_extension_array_equal(result, expected)
    # Test using series accessor
    expected = pd.Series(
        result, name=f"{data_for_encoding.variant.id}_{data_for_encoding.variant.alt}"
    )
    result_series = pd.Series(data_for_encoding).genotype.encode_additive()
    assert_series_equal(result_series, expected)


def test_encoding_dominant(data_for_encoding):
    # Test arrays directly
    expected = pd.array([0, 1, 1, None], dtype="UInt8")
    result = data_for_encoding.encode_dominant()
    assert_extension_array_equal(result, expected)
    # Test using series accessor
    expected = pd.Series(
        result, name=f"{data_for_encoding.variant.id}_{data_for_encoding.variant.alt}"
    )
    result_series = pd.Series(data_for_encoding).genotype.encode_dominant()
    assert_series_equal(result_series, expected)


def test_encoding_recessive(data_for_encoding):
    # Test arrays directly
    expected = pd.array([0, 0, 1, None], dtype="UInt8")
    result = data_for_encoding.encode_recessive()
    assert_extension_array_equal(result, expected)
    # Test using series accessor
    expected = pd.Series(
        result, name=f"{data_for_encoding.variant.id}_{data_for_encoding.variant.alt}"
    )
    result_series = pd.Series(data_for_encoding).genotype.encode_recessive()
    assert_series_equal(result_series, expected)


def test_encoding_codominant(data_for_encoding):
    # Test arrays directly
    expected = pd.Categorical(
        ["Ref", "Het", "Hom", None], categories=["Ref", "Het", "Hom"], ordered=True
    )
    result = data_for_encoding.encode_codominant()
    assert_extension_array_equal(result, expected)
    # Test using series accessor
    expected = pd.Series(
        result, name=f"{data_for_encoding.variant.id}_{data_for_encoding.variant.alt}"
    )
    result_series = pd.Series(data_for_encoding).genotype.encode_codominant()
    assert_series_equal(result_series, expected)


def test_var_info(genotypearray_df):
    print(genotypearray_df.genotyping.variant_info)
