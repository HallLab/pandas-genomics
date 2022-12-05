"""
Test GenotypeArray Accessors
"""
import numpy as np
import pandas as pd
import pytest
from pandas._testing import assert_frame_equal, assert_series_equal

from pandas_genomics import GenotypeArray, sim
from pandas_genomics.scalars import Region, Variant


def test_variant_score(data, data_for_encoding):
    assert pd.Series(data).genomics.variant.score == float(data.variant.score)
    assert pd.Series(data_for_encoding()).genomics.variant.score is None


def test_var_info(genotypearray_df):
    assert_series_equal(
        genotypearray_df.genomics.variant_info.iloc[0],
        genotypearray_df.iloc[:, 0].genomics.variant_info,
    )


def test_maf(data):
    assert pd.Series(data).genomics.maf == data.maf
    df = pd.DataFrame.from_dict(
        {n: pd.Series(data.copy()) for n in "ABC"}, orient="columns"
    )
    for colname in "ABC":
        df[colname].genomics.variant.id = colname
    df["D"] = np.ones(len(data))
    expected = pd.Series({"A": data.maf, "B": data.maf, "C": data.maf})
    assert_series_equal(df.genomics.maf, expected)


# def test_hwe(data):
#     assert pd.Series(data).genomics.hwe_pval == data.hwe_pval


@pytest.mark.parametrize(
    "filter_value, num_cols_left", [(None, 17), (0.05, 3), (0.10, 2)]
)
def test_filter_maf(genotypearray_df, filter_value, num_cols_left):
    if filter_value is None:
        result = genotypearray_df.genomics.filter_variants_maf()
    else:
        result = genotypearray_df.genomics.filter_variants_maf(filter_value)
    assert len(result.columns) == num_cols_left


@pytest.mark.parametrize(
    "filter_value, num_cols_left", [(None, 1), (0.05, 1), (1e-300, 2)]
)
def test_filter_hwe(ga_inhwe, ga_nothwe, filter_value, num_cols_left):
    data = pd.DataFrame({"yes": ga_inhwe, "no": ga_nothwe})
    data["num"] = [n for n in range(len(data))]
    if filter_value is None:
        result = data.genomics.filter_variants_hwe()
    else:
        result = data.genomics.filter_variants_hwe(filter_value)
    assert len(result.columns) == num_cols_left + 1


def test_region_series():
    var = Variant("chr1", position=999, ref="A", alt=["a"])
    s = pd.Series(
        GenotypeArray(
            [
                var.make_genotype_from_str("A/A"),
            ]
            * 10
        )
    )
    assert s.genomics.contained_by(Region("chr1", 900, 1000))
    assert not s.genomics.contained_by(Region("chr2", 900, 1000))
    assert not s.genomics.contained_by(Region("chr1", 900, 999))


def test_region_df():
    var1 = Variant("chr1", position=999, ref="A", alt=["a"])
    var2 = Variant("chr1", position=6789, ref="A", alt=["a"])
    var3 = Variant("chr2", position=999, ref="A", alt=["a"])
    var4 = Variant("chr3", position=25622, ref="A", alt=["a"])
    df = pd.DataFrame(
        {
            f"var{idx+1}": GenotypeArray(
                [
                    var.make_genotype_from_str("A/A"),
                ]
                * 10
            )
            for idx, var in enumerate([var1, var2, var3, var4])
        }
    )
    assert_frame_equal(
        df.genomics.in_regions(Region("chr1", 900, 1000)),
        df[
            [
                "var1",
            ]
        ],
    )
    assert_frame_equal(
        df.genomics.not_in_regions(Region("chr1", 900, 1000)),
        df[["var2", "var3", "var4"]],
    )
    assert_frame_equal(
        df.genomics.in_regions([Region("chr1", 900, 1000), Region("chr2", 900, 1000)]),
        df[["var1", "var3"]],
    )


def test_routine_case_control():
    # Additive Main Effect for SNP1 without interaction
    train = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=0.45,
        main1=1,
        main2=0,
        interaction=0,
        random_seed=2021,
    )
    train_cases = train.generate_case_control(n_cases=5000, n_controls=5000)
    edge = train_cases.genomics.calculate_edge_encoding_values(
        data=train_cases["Outcome"], outcome_variable="Outcome"
    )
    assert edge["Variant ID"][0] == "rs1"
