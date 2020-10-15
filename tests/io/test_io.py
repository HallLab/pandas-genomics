import pytest
import pandas_genomics

@pytest.mark.slow
def test_loaded_small(plink_small):
    """Validate the small dataset"""
    # TODO: Add more assertions
    assert plink_small.shape == (150, 3020)


@pytest.mark.slow
def test_loaded_medium(plink_medium):
    """Validate the medium dataset"""
    # TODO: Add more assertions
    assert plink_medium.shape == (600, 45100)


def test_loaded_small_20(plink_small_20):
    """Validate the small dataset limited to 20 variants"""
    # TODO: Add more assertions
    assert plink_small_20.shape == (150, 20)


def test_loaded_small_20_swap(plink_small_20_swap):
    """Validate the small dataset limited to 20 variants with swapped alleles"""
    # TODO: Add more assertions
    assert plink_small_20_swap.shape == (150, 20)


def test_swap(plink_small_20, plink_small_20_swap):
    """Validate that swapping alleles worked when loading"""
    normal = plink_small_20.iloc[:, 0]
    swapped = plink_small_20_swap.iloc[:, 0]
    assert normal.dtype != swapped.dtype
    swapped.genotype.set_reference(1)
    assert (normal == swapped).all()
