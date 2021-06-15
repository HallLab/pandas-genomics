"""
Test GenotypeArray Encoding methods and accessors to those functions
"""
import copy

import pandas as pd
import pytest
from pandas._testing import (
    assert_extension_array_equal,
    assert_series_equal,
    assert_frame_equal,
)


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
    result_series = pd.Series(data_for_encoding).genomics.encode_additive()
    assert_series_equal(result_series, expected)
    # Test using DataFrame accessor
    df = pd.DataFrame.from_dict({n: data_for_encoding for n in "ABC"}, orient="columns")
    expected = pd.concat([result_series] * 3, axis=1)
    result_df = df.genomics.encode_additive()
    assert_frame_equal(result_df, expected)


def test_encoding_dominant(data_for_encoding):
    # Test arrays directly
    expected = pd.array([0, 1, 1, None], dtype="UInt8")
    result = data_for_encoding.encode_dominant()
    assert_extension_array_equal(result, expected)
    # Test using series accessor
    expected = pd.Series(
        result, name=f"{data_for_encoding.variant.id}_{data_for_encoding.variant.alt}"
    )
    result_series = pd.Series(data_for_encoding).genomics.encode_dominant()
    assert_series_equal(result_series, expected)
    # Test using DataFrame accessor
    df = pd.DataFrame.from_dict({n: data_for_encoding for n in "ABC"}, orient="columns")
    expected = pd.concat([result_series] * 3, axis=1)
    result_df = df.genomics.encode_dominant()
    assert_frame_equal(result_df, expected)


def test_encoding_recessive(data_for_encoding):
    # Test arrays directly
    expected = pd.array([0, 0, 1, None], dtype="UInt8")
    result = data_for_encoding.encode_recessive()
    assert_extension_array_equal(result, expected)
    # Test using series accessor
    expected = pd.Series(
        result, name=f"{data_for_encoding.variant.id}_{data_for_encoding.variant.alt}"
    )
    result_series = pd.Series(data_for_encoding).genomics.encode_recessive()
    assert_series_equal(result_series, expected)
    # Test using DataFrame accessor
    df = pd.DataFrame.from_dict({n: data_for_encoding for n in "ABC"}, orient="columns")
    expected = pd.concat([result_series] * 3, axis=1)
    result_df = df.genomics.encode_recessive()
    assert_frame_equal(result_df, expected)


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
    result_series = pd.Series(data_for_encoding).genomics.encode_codominant()
    assert_series_equal(result_series, expected)
    # Test using DataFrame accessor
    df = pd.DataFrame.from_dict({n: data_for_encoding for n in "ABC"}, orient="columns")
    expected = pd.concat([result_series] * 3, axis=1)
    result_df = df.genomics.encode_codominant()
    assert_frame_equal(result_df, expected)


@pytest.mark.parametrize(
    "alpha_value,ref_allele,alt_allele,minor_allele_freq,expected",
    [
        (0.25, "A", "T", 0.45, pd.array([0.0, 0.25, 1.0, None], dtype="Float64")),
        (0.25, "T", "A", 0.45, pd.array([1.0, 0.25, 0.0, None], dtype="Float64")),
        (1.0, "T", "A", 0.45, pd.array([1.0, 1.0, 0.0, None], dtype="Float64")),
        pytest.param(0.25, "A", "C", 0.45, None,
                     marks=pytest.mark.xfail(raises=ValueError, strict=True, reason="Wrong Allele"))
    ],
)
def test_encoding_weighted(data_for_encoding, alpha_value, ref_allele, alt_allele, minor_allele_freq, expected):
    result = data_for_encoding.encode_weighted(alpha_value, ref_allele, alt_allele, minor_allele_freq)
    assert_extension_array_equal(expected, result)


def test_encoding_weighted_df(data_for_encoding):
    df = pd.DataFrame({name: copy.deepcopy(data_for_encoding) for name in ["Var1", "Var2", "Var3"]})
    for colname, col in df.iteritems():
        col.array.variant.id = colname
    encoding_info = pd.DataFrame({"Variant ID": ["Var1", "Var2", "Var3"],
                                  "Alpha Value": [0.70, 1.14, 0.21],
                                  "Ref Allele": ["T", "T", "T"],
                                  "Alt Allele": ["A", "A", "A"],
                                  "Minor Allele Frequency": [0.1, 0.1, 0.1]})
    result = df.genomics.encode_weighted(encoding_info)
    assert_extension_array_equal(expected, result)
