from pathlib import Path

import pytest

from pandas_genomics import io

DATA_DIR = Path(__file__).parent.parent / "data" / "plink"

# TODO: Check genotypes in plink and double-check that the results are correct


def test_loaded_small():
    """Validate the small dataset"""
    # was 5.76 seconds
    bed_file = DATA_DIR / "plink_test_small.bed"
    result = io.from_plink(bed_file)
    assert result.shape == (150, 3020)


@pytest.mark.slow
def test_loaded_medium():
    """Validate the medium dataset"""
    bed_file = DATA_DIR / "plink_test_medium.bed"
    result = io.from_plink(bed_file)
    assert result.shape == (600, 45100)


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
    swapped.genomics.set_reference(1)
    assert normal.genomics.variant == swapped.genomics.variant
    assert (normal == swapped).all()
