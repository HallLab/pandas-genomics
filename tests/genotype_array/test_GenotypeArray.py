"""
Test GenotypeArray methods
"""
import pytest


# Custom Tests
@pytest.mark.xfail(raises=ValueError)
def test_encoding_additive(data):
    result = data.encode_additive()
