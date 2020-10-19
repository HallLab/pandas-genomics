import sys
from pathlib import Path

import pytest

from pandas_genomics import io

DATA_DIR = Path(__file__).parent.parent / "data"


@pytest.mark.slow
@pytest.fixture
def plink_small():
    bed_file = DATA_DIR / "plink" / "plink_test_small.bed"
    result = io.from_plink(bed_file)
    return result


@pytest.mark.slow
@pytest.fixture
def plink_medium():
    bed_file = DATA_DIR / "plink" / "plink_test_medium.bed"
    result = io.from_plink(bed_file)
    return result


@pytest.fixture
def plink_small_20():
    bed_file = DATA_DIR / "plink" / "plink_test_small.bed"
    result = io.from_plink(bed_file, max_variants=20)
    return result


@pytest.fixture
def plink_small_20_swap():
    bed_file = DATA_DIR / "plink" / "plink_test_small.bed"
    result = io.from_plink(bed_file, max_variants=20, swap_alleles=True)
    return result


@pytest.fixture
def vcf_test():
    # Skip on Windows
    if sys.platform.startswith("win"):
        pytest.skip("VCF IO requires HTSLIB, which isn't easy to install on Windows")
    # Run otherwise
    import cyvcf2

    VCF_DIR = Path(cyvcf2.__file__).parent / "tests"
    vcf_filename = VCF_DIR / "test.vcf.gz"
    result = io.from_vcf(vcf_filename)
    return result
