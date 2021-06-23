"""
Test GenotypeDtype
"""
import pandas as pd
import pytest
from pandas._testing import assert_series_equal, assert_extension_array_equal

from pandas_genomics.arrays import GenotypeDtype
from pandas_genomics.scalars import Variant

TEST_VAR = Variant(
    chromosome="12", position=112161652, id="rs12462", ref="T", alt=["C"]
)


@pytest.mark.parametrize(
    "input_str,variant",
    [
        (
            "genotype(2n)[12; 112161652; rs12462; T; C]",
            Variant(
                chromosome="12", position=112161652, id="rs12462", ref="T", alt=["C"]
            ),
        ),
        (
            "genotype(3n)[12; 112161652; rs12462; T; C]",
            Variant(
                chromosome="12",
                position=112161652,
                id="rs12462",
                ref="T",
                alt=["C"],
                ploidy=3,
            ),
        ),
        (
            "genotype(3n)[12; 112161652; rs12462; T; C]Q30",
            Variant(
                chromosome="12",
                position=112161652,
                id="rs12462",
                ref="T",
                alt=["C"],
                ploidy=3,
                score=30,
            ),
        ),
        pytest.param(
            "genotype[12; 112161652; rs12462; T; C]",
            None,
            marks=pytest.mark.xfail(raises=TypeError, strict=True),
        ),
        pytest.param(
            "genotype(2n)[12; 112161652; T; C]",
            None,
            marks=pytest.mark.xfail(raises=TypeError, strict=True),
        ),
        pytest.param(
            "genotype(2n)[12; 112161652; T; C]q35",
            None,
            marks=pytest.mark.xfail(raises=TypeError, strict=True),
        ),
    ],
)
def test_from_str(input_str, variant):
    """Test creating GenotypeDtype from str"""
    gtdtype = GenotypeDtype.construct_from_string(input_str)
    assert gtdtype.variant == variant
