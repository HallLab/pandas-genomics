import random

from pandas.tests.extension import base
import pytest

from pandas_genomics.arrays import GenotypeDtype, GenotypeArray
from pandas_genomics.scalars import Variant

random.seed(1855)


# Implement the required fixtures
@pytest.fixture
def dtype():
    variant = Variant(variant_id='rs12345', chromosome='chr1', coordinate=123456, alleles=['A', 'T', 'G'])
    return GenotypeDtype(variant=variant)


@pytest.fixture
def data():
    """Length-100 array for this type.
    * data[0] and data[1] should both be non missing
    * data[0] and data[1] should not be equal
    """
    alleles = ['A', 'T', 'G']
    variant = Variant(variant_id='rs12345', chromosome='chr1', coordinate=123456, alleles=alleles)
    genotypes = [variant.make_genotype('A', 'T'), variant.make_genotype('T', 'T')]
    for i in range(98):
        genotypes.append(variant.make_genotype(random.choice(alleles), random.choice(alleles)))
    return GenotypeArray(values=genotypes)


@pytest.fixture
def data_for_twos():
    """Length-100 array in which all the elements are two."""
    # Not applicable
    raise NotImplementedError


@pytest.fixture
def data_missing():
    """Length-2 array with [NA, Valid]"""
    alleles = ['A', 'T', 'G']
    variant = Variant(variant_id='rs12345', chromosome='chr1', coordinate=123456, alleles=alleles)
    genotypes = [variant.make_genotype(), variant.make_genotype('T', 'T')]
    return GenotypeArray(values=genotypes)


@pytest.fixture
def data_for_sorting():
    """Length-3 array with a known sort order.
    This should be three items [B, C, A] with
    A < B < C
    """
    alleles = ['A', 'T', 'G']
    variant = Variant(variant_id='rs12345', chromosome='chr1', coordinate=123456, alleles=alleles)
    a = variant.make_genotype('A', 'A')
    b = variant.make_genotype('A', 'T')
    c = variant.make_genotype('T', 'T')
    return GenotypeArray(values=[b, c, a])


@pytest.fixture
def data_missing_for_sorting():
    """Length-3 array with a known sort order.
    This should be three items [B, NA, A] with
    A < B and NA missing.
    """
    alleles = ['A', 'T', 'G']
    variant = Variant(variant_id='rs12345', chromosome='chr1', coordinate=123456, alleles=alleles)
    a = variant.make_genotype('A', 'A')
    b = variant.make_genotype('A', 'T')
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
    variant = Variant(variant_id='rs12345', chromosome='chr1', coordinate=123456, alleles=['A', 'T', 'G'])
    return variant.make_genotype()


@pytest.fixture
def data_for_grouping():
    """Data for factorization, grouping, and unique tests.
    Expected to be like [B, B, NA, NA, A, A, B, C]
    Where A < B < C and NA is missing
    """
    alleles = ['A', 'T', 'G']
    variant = Variant(variant_id='rs12345', chromosome='chr1', coordinate=123456, alleles=alleles)
    a = variant.make_genotype('A', 'A')
    b = variant.make_genotype('A', 'T')
    c = variant.make_genotype('T', 'T')
    na = variant.make_genotype()
    return GenotypeArray([b, b, na, na, a, a, b, c])


# Run the predefined tests
class TestCasting(base.BaseCastingTests):
    pass


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestDtype(base.BaseDtypeTests):
    pass


class TestGetItem(base.BaseGetitemTests):
    pass


class TestGroupBy(base.BaseGroupbyTests):
    pass


class TestInterface(base.BaseInterfaceTests):
    pass


class TestParsing(base.BaseParsingTests):
    pass


class TestMethods(base.BaseMethodsTests):

    def test_combine_add(self, data_repeated):
        """Addition of Genotypes isn't valid"""
        pass

    def test_searchsorted(self, data_for_sorting, as_series):
        # TODO: Can't pass until it's possible to define dtype as scalar (See Pandas GH #33825)
        pass

    def test_where_series(self, data, na_value, as_frame):
        # TODO: Can't pass until it's possible to define dtype as scalar (See Pandas GH #33825)
        pass


class TestMissing(base.BaseMissingTests):
    pass


# Skip ArithmeticOps since they aren't valid
# class TestArithmeticOps(base.BaseArithmeticOpsTests):
#     pass


class TestComparisonOps(base.BaseComparisonOpsTests):
    pass


class TestOpsUtil(base.BaseOpsUtil):
    pass

# No way to invert a genotype
# class TestUnaryOps(base.BaseUnaryOpsTests):
#     pass


class TestPrinting(base.BasePrintingTests):
    pass


# No boolean equivalent for genotypes
# class TestBooleanReduce(base.BaseBooleanReduceTests):
#     pass


class TestNoReduce(base.BaseNoReduceTests):
    pass


# No numeric equivalent for genotypes
# class TestNumericReduce(base.BaseNumericReduceTests):
#     pass


class TestReshaping(base.BaseReshapingTests):
    pass


class TestSetitems(base.BaseSetitemTests):
    pass
