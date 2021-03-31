"""
Test various calculations performed by GenotypeArray
"""
import numpy as np

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


def test_maf():
    var = Variant("chr1", ref="A", alt=["T", "C"])
    ga_zero = GenotypeArray(
        [
            var.make_genotype_from_str("A/A"),
            var.make_genotype_from_str("A/A"),
            var.make_genotype_from_str("A/A"),
        ]
    )
    assert ga_zero.maf == 0.0
    ga_50 = GenotypeArray(
        [
            var.make_genotype_from_str("A/A"),
            var.make_genotype_from_str("T/T"),
            var.make_genotype_from_str("T/A"),
        ]
    )
    assert ga_50.maf == 0.50
    ga_2nd = GenotypeArray(
        [
            var.make_genotype_from_str("A/C"),
            var.make_genotype_from_str("C/C"),
            var.make_genotype_from_str("T/T"),
        ]
    )
    assert ga_2nd.maf == 0.5
    missing = GenotypeArray([var.make_genotype()] * 3)
    assert missing.maf is np.nan
