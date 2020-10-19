"""
Run ExtensionArray tests from Pandas on the GenotypeArray class
"""

from pandas.tests.extension import base


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
    def test_unstack(self, data, index, obj):
        # Fixed for Pandas 1.2: https://github.com/pandas-dev/pandas/issues/36986
        pass


class TestSetitems(base.BaseSetitemTests):
    pass
