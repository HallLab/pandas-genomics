"""
Test GenotypeDtype
"""
import pandas as pd
import pytest
from pandas._testing import assert_series_equal, assert_extension_array_equal

from pandas_genomics.arrays import GenotypeDtype
from pandas_genomics.scalars import Variant

TEST_VAR = Variant(chromosome="12", position=112161652, id="rs12462", ref="T", alt=["C"])


@pytest.mark.parametrize(
    "input_str,variant,ploidy",
    [("genotype(2)[12; 112161652; rs12462; T; C]",
      Variant(chromosome="12", position=112161652, id="rs12462", ref="T", alt=["C"]), 2),
     ("genotype(3)[12; 112161652; rs12462; T; C]",
      Variant(chromosome="12", position=112161652, id="rs12462", ref="T", alt=["C"]), 3),
     pytest.param("genotype[12; 112161652; rs12462; T; C]", None, None, marks=pytest.mark.xfail(raises=TypeError)),
     pytest.param("genotype(2)[12; 112161652; T; C]", None, None, marks=pytest.mark.xfail(raises=TypeError)),
     pytest.param("genotype(2)[12; 112161652; rs12462; T; C]",
                  Variant(chromosome="12", position=112161652, id="rs12462", ref="T", alt=["C"]),
                  3,
                  marks=pytest.mark.xfail(raises=AssertionError))
     ])
def test_from_str(input_str, variant, ploidy):
    """Test creating GenotypeDtype from str"""
    gtdtype = GenotypeDtype.construct_from_string(input_str)
    assert gtdtype.variant == variant
    assert gtdtype.ploidy == ploidy
