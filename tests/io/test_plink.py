from pathlib import Path

import pytest
from numpy.testing import assert_array_equal
from pandas._testing import assert_frame_equal

from pandas_genomics import io, sim

DATA_DIR = Path(__file__).parent.parent / "data" / "plink"

# TODO: Check genotypes in plink and double-check that the results are correct


def test_small():
    """Validate the small dataset"""
    input = DATA_DIR / "plink_test_small"
    result = io.from_plink(input, categorical_phenotype=True)
    assert result.shape == (150, 3020)


def test_round_trip_real(tmp_path):
    """Load real data, save it, and load it again"""
    d = tmp_path / "test"
    d.mkdir()
    output = str(d / "test")
    # Load data
    input = DATA_DIR / "plink_test_small"
    loaded = io.from_plink(str(input), categorical_phenotype=True, max_variants=100)
    # Save data
    io.to_plink(loaded, output)
    # Reload data
    reloaded = io.from_plink(str(output), categorical_phenotype=True)
    # Compare
    assert_frame_equal(
        loaded.reset_index(), reloaded.reset_index(), check_categorical=False
    )


def test_round_trip_sim(tmp_path):
    """Simulate data, save it, and load it again"""
    d = tmp_path / "test"
    d.mkdir()
    output = str(d / "test")
    data = sim.BAMS().generate_case_control()
    original = data.copy()
    io.to_plink(
        data,
        output,
        phenotype_name="Outcome",
        phenotype_case="Case",
        phenotype_control="Control",
    )
    # Load data and reset index to extract phenotype and get original data format back
    loaded_data = (
        io.from_plink(output, categorical_phenotype=True)
        .reset_index(level=-1)
        .reset_index(drop=True)
    )
    loaded_data.columns = data.columns  # Correct column names
    assert_frame_equal(
        data, loaded_data, check_categorical=False
    )  # Categorical order may be different

    # Ensure there were no side effects
    assert_array_equal(
        original["SNP1"].array.allele_idxs, data["SNP1"].array.allele_idxs
    )


@pytest.mark.slow
def test_loaded_medium():
    """Validate the medium dataset"""
    input = DATA_DIR / "plink_test_medium"
    result = io.from_plink(input)
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
