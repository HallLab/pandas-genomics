from pathlib import Path
import pytest

from pandas_genomics import io

data_dir = Path(__file__).parent.parent / "data"


@pytest.mark.slow
@pytest.fixture
def plink_small():
    bed_file = data_dir / "plink" / "plink_test_small.bed"
    result = io.from_plink(bed_file)
    return result


@pytest.mark.slow
@pytest.fixture
def plink_medium():
    bed_file = data_dir / "plink" / "plink_test_medium.bed"
    result = io.from_plink(bed_file)
    return result


@pytest.fixture
def plink_small_20():
    bed_file = data_dir / "plink" / "plink_test_small.bed"
    result = io.from_plink(bed_file, max_variants=20)
    return result


@pytest.fixture
def plink_small_20_swap():
    bed_file = data_dir / "plink" / "plink_test_small.bed"
    result = io.from_plink(bed_file, max_variants=20, swap_alleles=True)
    return result
