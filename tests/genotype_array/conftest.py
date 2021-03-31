import random

import pytest

import pandas as pd
from pandas_genomics.arrays import GenotypeDtype, GenotypeArray
from pandas_genomics.scalars import Variant

random.seed(1855)


# Below fixtures are copied from pandas.conftest
# They could be imported, but that would require having hypothesis as a dependency
@pytest.fixture(params=[None, lambda x: x])
def sort_by_key(request):
    """
    Simple fixture for testing keys in sorting methods.
    Tests None (no key) and the identity key.
    """
    return request.param


@pytest.fixture(params=["__eq__", "__ne__", "__le__", "__lt__", "__ge__", "__gt__"])
def all_compare_operators(request):
    """
    Fixture for dunder names for common compare operations
    * >=
    * >
    * ==
    * !=
    * <
    * <=
    """
    return request.param


_all_numeric_reductions = [
    "sum",
    "max",
    "min",
    "mean",
    "prod",
    "std",
    "var",
    "median",
    "kurt",
    "skew",
]


@pytest.fixture(params=_all_numeric_reductions)
def all_numeric_reductions(request):
    """
    Fixture for numeric reduction names.
    """
    return request.param


_all_boolean_reductions = ["all", "any"]


@pytest.fixture(params=_all_boolean_reductions)
def all_boolean_reductions(request):
    """
    Fixture for boolean reduction names.
    """
    return request.param


# Implement the required fixtures
@pytest.fixture
def dtype():
    variant = Variant(
        chromosome="chr1",
        position=123456,
        id="rs12345",
        ref="A",
        alt=["T", "G"],
        score=30,
    )
    return GenotypeDtype(variant=variant)


@pytest.fixture
def data():
    """Length-100 array for this type.
    * data[0] and data[1] should both be non missing
    * data[0] and data[1] should not be equal
    """
    alleles = ["A", "T", "G"]
    variant = Variant(
        chromosome="chr1",
        position=123456,
        id="rs12345",
        ref="A",
        alt=["T", "G"],
        score=30,
    )
    genotypes = [variant.make_genotype("A", "T"), variant.make_genotype("T", "T")]
    for i in range(98):
        genotypes.append(
            variant.make_genotype(random.choice(alleles), random.choice(alleles))
        )
    return GenotypeArray(values=genotypes)


@pytest.fixture
def data_for_twos():
    """Length-100 array in which all the elements are two."""
    # Not applicable
    raise NotImplementedError


@pytest.fixture
def data_missing():
    """Length-2 array with [NA, Valid]"""
    variant = Variant(
        chromosome="chr1", position=123456, id="rs12345", ref="A", alt=["T", "G"]
    )
    genotypes = [variant.make_genotype(), variant.make_genotype("T", "T")]
    return GenotypeArray(values=genotypes)


@pytest.fixture
def data_for_sorting():
    """Length-3 array with a known sort order.
    This should be three items [B, C, A] with
    A < B < C
    """
    variant = Variant(
        chromosome="chr1",
        position=123456,
        id="rs12345",
        ref="A",
        alt=["T", "G"],
        score=30,
    )
    a = variant.make_genotype("A", "A")
    b = variant.make_genotype("A", "T")
    c = variant.make_genotype("T", "T")
    return GenotypeArray(values=[b, c, a])


@pytest.fixture
def data_missing_for_sorting():
    """Length-3 array with a known sort order.
    This should be three items [B, NA, A] with
    A < B and NA missing.
    """
    variant = Variant(
        chromosome="chr1", position=123456, id="rs12345", ref="A", alt=["T", "G"]
    )
    a = variant.make_genotype("A", "A")
    b = variant.make_genotype("A", "T")
    na = variant.make_genotype()
    return GenotypeArray(values=[b, na, a])


@pytest.fixture
def na_cmp():
    """Binary operator for comparing NA values.
    Should return a function of two arguments that returns
    True if both arguments are (scalar) NA for your type.
    By default, uses ``operator.is_``
    """
    return lambda gt1, gt2: gt1.is_missing() and gt2.is_missing()


@pytest.fixture
def na_value():
    """The scalar missing value for this type. Default 'None'"""
    variant = Variant(
        chromosome="chr1",
        position=123456,
        id="rs12345",
        ref="A",
        alt=["T", "G"],
        score=30,
    )
    return variant.make_genotype()


@pytest.fixture
def data_for_grouping():
    """Data for factorization, grouping, and unique tests.
    Expected to be like [B, B, NA, NA, A, A, B, C]
    Where A < B < C and NA is missing
    """
    variant = Variant(
        chromosome="chr1", position=123456, id="rs12345", ref="A", alt=["T", "G"]
    )
    a = variant.make_genotype("A", "A")
    b = variant.make_genotype("A", "T")
    c = variant.make_genotype("T", "T")
    na = variant.make_genotype()
    return GenotypeArray([b, b, na, na, a, a, b, c])


@pytest.fixture
def data_for_encoding():
    """Data for encoding tests.
    Contains one alt allele.
    Variants are Homozygouse Ref, Heterozygous, Homozygous Alt, and Missing
    """
    variant = Variant(
        chromosome="chr1", position=123456, id="rs12345", ref="A", alt=["T"]
    )
    a = variant.make_genotype("A", "A")
    b = variant.make_genotype("A", "T")
    c = variant.make_genotype("T", "T")
    na = variant.make_genotype()
    return GenotypeArray([a, b, c, na])

@pytest.fixture
def genotypearray_df():
    variants = [
        Variant(chromosome="chr1", position=123456, id="rs12345", ref="A", alt=["T"], score=30),
        Variant(chromosome="chr1", position=123456, id="rs12345", ref="T", alt=["A"], score=30),
        Variant(chromosome="chr2", position=130546, id="rs12345", ref="C", alt=["G", "T"], score=30),
        Variant(chromosome="chr2", position=131111, id="rs12345", ref="A", alt=["T"], score=30),
        Variant(chromosome="chrX", position=34245, id="rs12345", ref="A", alt=["T"], score=30),
    ]
    return pd.DataFrame({name:
        GenotypeArray([var.make_genotype(var.ref, var.alt[0]),
                       var.make_genotype(var.ref, var.ref),
                       var.make_genotype(var.alt[0], var.alt[0]),
                       var.make_genotype(var.ref, var.ref),
                       var.make_genotype(var.ref, var.alt[0])])
                      for var, name in zip(variants, "ABCDE")})
