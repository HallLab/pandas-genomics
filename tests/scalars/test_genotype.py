import pytest

from pandas_genomics.scalars import Genotype, Variant


@pytest.fixture()
def var_complete():
    return Variant("12", 12345678, "complete", ref="A", alt=["C", "G", "T"])


@pytest.fixture()
def var_two_allele():
    return Variant("12", 12345678, "complete", ref="A", alt=["C"])


@pytest.fixture()
def var_min():
    return Variant("13", 12345678, "min", ref="A", alt=[])


def test_make_genotype(var_complete):
    gt = var_complete.make_genotype("A", "C")
    expected = Genotype(variant=var_complete, allele_idxs=[0, 1])
    assert gt == expected
    # Triploid and sorting alleles correctly
    gt = var_complete.make_genotype("A", "C", "A")
    expected = Genotype(variant=var_complete, allele_idxs=[0, 0, 1])
    assert gt == expected


def test_missing(var_complete):
    gt = var_complete.make_genotype()
    assert gt.is_missing()


def test_make_genotype_from_str(var_complete):
    # Haploid
    gt = var_complete.make_genotype_from_str("A")
    expected = Genotype(variant=var_complete, allele_idxs=[0])
    assert gt == expected
    # Diploid
    gt = var_complete.make_genotype_from_str("A/C")
    expected = Genotype(variant=var_complete, allele_idxs=[0, 1])
    assert gt == expected
    # Triploid and sorting
    gt = var_complete.make_genotype_from_str("A/C/A")
    expected = Genotype(variant=var_complete, allele_idxs=[0, 0, 1])
    assert gt == expected


def test_make_genotype_from_plink_bits(var_two_allele):
    gt = var_two_allele.make_genotype_from_plink_bits("00")
    expected = Genotype(variant=var_two_allele, allele_idxs=[0, 0])
    assert gt == expected
    gt = var_two_allele.make_genotype_from_plink_bits("01")
    expected = Genotype(variant=var_two_allele)
    assert gt == expected
    gt = var_two_allele.make_genotype_from_plink_bits("10")
    expected = Genotype(variant=var_two_allele, allele_idxs=[0, 1])
    assert gt == expected
    gt = var_two_allele.make_genotype_from_plink_bits("11")
    expected = Genotype(variant=var_two_allele, allele_idxs=[1, 1])
    assert gt == expected
