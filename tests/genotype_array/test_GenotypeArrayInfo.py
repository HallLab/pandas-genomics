"""
Test various calculations performed by GenotypeArray
"""
import numpy as np
import pandas as pd

from pandas_genomics.arrays import GenotypeArray
from pandas_genomics.scalars import Variant


def test_is_missing():
    var = Variant("chr1", ref="A", alt=["T", "C"])
    ga_fff = GenotypeArray(
        [
            var.make_genotype_from_str("A/A"),
            var.make_genotype_from_str("A/A"),
            var.make_genotype_from_str("A/A"),
        ]
    )
    assert (ga_fff.is_missing == np.array([False, False, False])).all()
    ga_ftf = GenotypeArray(
        [
            var.make_genotype_from_str("A/A"),
            var.make_genotype(),
            var.make_genotype_from_str("A/A"),
        ]
    )
    assert (ga_ftf.is_missing == np.array([False, True, False])).all()
    assert (ga_ftf.is_missing == (ga_ftf == "./.")).all()


def test_equality(ga_AA_Aa_aa_BB_Bb_bb):
    assert (
        (ga_AA_Aa_aa_BB_Bb_bb == "a/a")
        == np.array([False, False, True, False, False, False])
    ).all()
    assert (
        (pd.Series(ga_AA_Aa_aa_BB_Bb_bb) == "a/a")
        == np.array([False, False, True, False, False, False])
    ).all()


def test_is_funcs(ga_AA_Aa_aa_BB_Bb_bb):
    assert (
        ga_AA_Aa_aa_BB_Bb_bb.is_homozygous
        == np.array([True, False, True, True, False, True])
    ).all()
    assert (
        ga_AA_Aa_aa_BB_Bb_bb.is_heterozygous
        == np.array([False, True, False, False, True, False])
    ).all()
    assert (
        ga_AA_Aa_aa_BB_Bb_bb.is_homozygous_ref
        == np.array([True, False, False, False, False, False])
    ).all()
    assert (
        ga_AA_Aa_aa_BB_Bb_bb.is_homozygous_alt
        == np.array([False, False, True, True, False, True])
    ).all()


def test_maf():
    # Zero
    var = Variant("chr1", ref="A", alt=["T", "C"])
    ga_zero = GenotypeArray(
        [
            var.make_genotype_from_str("A/A"),
            var.make_genotype_from_str("A/A"),
            var.make_genotype_from_str("A/A"),
        ]
    )
    assert ga_zero.maf == 0.0

    # Only Missing
    missing = GenotypeArray([var.make_genotype()] * 3)
    assert missing.maf is np.nan

    # 50%
    ga_50 = GenotypeArray(
        [
            var.make_genotype_from_str("A/A"),
            var.make_genotype_from_str("T/T"),
            var.make_genotype_from_str("T/A"),
        ]
    )
    assert ga_50.maf == 0.50

    # 2nd of 3 alleles
    ga_2nd = GenotypeArray(
        [
            var.make_genotype_from_str("A/C"),
            var.make_genotype_from_str("C/C"),
            var.make_genotype_from_str("T/T"),
        ]
    )
    assert ga_2nd.maf == 0.5

    # Triploid
    var = Variant("chr1", ref="A", alt=["T", "C"], ploidy=3)
    ga_33 = GenotypeArray(
        [
            var.make_genotype_from_str("A/A/T"),
            var.make_genotype_from_str("A/T/C"),
            var.make_genotype_from_str("A/A/T"),
        ]
    )
    assert ga_33.maf == 1 / 3


def test_HWE(ga_inhwe, ga_nothwe):
    var = Variant("chr1", ref="A", alt=["a"])
    # One var, can't calculate
    ga_onevar = GenotypeArray(
        [
            var.make_genotype_from_str("A/A"),
        ]
    )
    assert ga_onevar.hwe_pval is np.nan
    assert ga_inhwe.hwe_pval == 1.0
    assert ga_nothwe.hwe_pval < 1e-20

    # NaN for non-diploid
    var = Variant("chr1", ref="A", alt=["B", "C"], ploidy=3)
    ga_triploid = GenotypeArray(
        [
            var.make_genotype_from_str("A/A/A"),
        ]
        * 50
        + [
            var.make_genotype_from_str("A/A/B"),
        ]
        * 50,
    )
    assert ga_triploid.hwe_pval is np.nan
