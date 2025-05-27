import sys
from pathlib import Path

import pytest

from pandas_genomics import io

DATA_DIR = Path(__file__).parent.parent / "data"


@pytest.fixture
def plink_small_20():
    input = DATA_DIR / "plink" / "plink_test_small"
    result = io.from_plink(input, max_variants=20)
    return result


@pytest.fixture
def plink_small_20_swap():
    input = DATA_DIR / "plink" / "plink_test_small"
    result = io.from_plink(input, max_variants=20, swap_alleles=True)
    return result


@pytest.fixture
def vcf_test():
    # Skip on Windows
    if sys.platform.startswith("win"):
        pytest.skip("VCF IO requires HTSLIB, which isn't easy to install on Windows")
    # Run otherwise
    # Andre: Update to Python >= 3.10
    # import cyvcf2
    # VCF_DIR = Path(cyvcf2.__file__).parent / "tests"
    # vcf_filename = VCF_DIR / "test.vcf.gz"
    vcf_filename = DATA_DIR / "vcf" / "test.vcf"
    result = io.from_vcf(vcf_filename)
    return result
