import pytest

from pandas_genomics.scalars import Genotype, Variant

@pytest.fixture()
def var_complete():
    return Variant("12", 12345678, "complete", ref="A", alt=["C", "G", "T"])

@pytest.fixture()
def var_min():
    return Variant("13", 12345678, "min", ref="A", alt=[])


def test_make_genotype(var_complete):
    gt1 = var_complete.make_genotype("A", "C")
    expected = Genotype(variant=var_complete, allele_idxs=[0, 1])
    assert gt1 == expected


