"""
Test GenotypeArray methods
"""
import pandas as pd
import pytest


@pytest.mark.xfail(raises=ValueError)
def test_encoding_extra_alt(data):
    """Fail encoding when there are multiple alt alleles"""
    data.encode_additive()


def test_encoding_additive(data_for_encoding):
    expected = pd.array([0, 1, 2, None], dtype='UInt8')
    result = data_for_encoding.encode_additive()
    assert (result == expected).all()
    # Test using series accessor
    result_series = pd.Series(data_for_encoding).genotype.encode_additive()
    assert (pd.Series(result) == result_series).all()


def test_encoding_dominant(data_for_encoding):
    expected = pd.array([0, 1, 1, None], dtype='UInt8')
    result = data_for_encoding.encode_dominant()
    assert (result == expected).all()
    # Test using series accessor
    result_series = pd.Series(data_for_encoding).genotype.encode_dominant()
    assert (pd.Series(result) == result_series).all()


def test_encoding_recessive(data_for_encoding):
    expected = pd.array([0, 0, 1, None], dtype='UInt8')
    result = data_for_encoding.encode_recessive()
    assert (result == expected).all()
    # Test using series accessor
    result_series = pd.Series(data_for_encoding).genotype.encode_recessive()
    assert (pd.Series(result) == result_series).all()
