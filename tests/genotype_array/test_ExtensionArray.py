"""
Run ExtensionArray tests from Pandas on the GenotypeArray class
"""
import pytest
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
    @pytest.mark.xfail(reason="GH#39098: Converts agg result to object")
    def test_groupby_agg_extension(self, data_for_grouping):
        super().test_groupby_agg_extension(data_for_grouping)


class TestInterface(base.BaseInterfaceTests):
    @pytest.mark.xfail(reason="'contains' method not implemented")
    def test_contains(self, data, data_missing):
        super().test_contains(data, data_missing)


class TestParsing(base.BaseParsingTests):
    pass


class TestMethods(base.BaseMethodsTests):
    @pytest.mark.xfail(reason="Addition of Genotypes isn't valid")
    def test_combine_add(self, data_repeated):
        super().test_combine_add(data_repeated)

    @pytest.mark.xfail(
        reason="Can't pass until it's possible to define dtype as scalar (See Pandas GH #33825)"
    )
    def test_searchsorted(self, data_for_sorting, as_series):
        super().test_searchsorted(data_for_sorting, as_series)

    @pytest.mark.xfail(
        reason="Can't pass until it's possible to define dtype as scalar (See Pandas GH #33825)"
    )
    def test_where_series(self, data, na_value, as_frame):
        super().test_where_series(data, na_value, as_frame)


class TestMissing(base.BaseMissingTests):
    pass


# Skip ArithmeticOps since they aren't valid
# class TestArithmeticOps(base.BaseArithmeticOpsTests):
#     pass


# class TestComparisonOps(base.BaseComparisonOpsTests):
#     pass


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
