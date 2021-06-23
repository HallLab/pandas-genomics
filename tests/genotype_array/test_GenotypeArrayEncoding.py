"""
Test GenotypeArray Encoding methods and accessors to those functions
"""

import pandas as pd
import pytest
from pandas._testing import (
    assert_extension_array_equal,
    assert_series_equal,
    assert_frame_equal,
)

import numpy as np
from pandas_genomics import GenotypeArray
from pandas_genomics.scalars import Variant
from pandas_genomics.sim import BAMS, SNPEffectEncodings


@pytest.fixture
def encoding_df():
    """
    5 variants, 5 genotypes each:
      Homozygous Ref
      Heterozygous
      Homozygous Alt
      Missing one allele
      Missing both alleles
    """
    data = dict()
    for idx, base in enumerate("ABCDE"):
        var = Variant(
            chromosome="chr1",
            position=idx + 1,
            id=f"rs{idx+1}",
            ref=base,
            alt=[base.lower()],
        )
        data[f"var{idx}"] = GenotypeArray(
            [
                var.make_genotype(base, base),
                var.make_genotype(base, base.lower()),
                var.make_genotype(base.lower(), base.lower()),
                var.make_genotype(base),
                var.make_genotype(),
            ]
        )
    return pd.DataFrame(data)


def test_encoding_additive(data_for_encoding):
    # Test arrays directly
    expected = pd.array([0, 1, 2, 2, None], dtype="UInt8")
    result = data_for_encoding().encode_additive()
    assert_extension_array_equal(result, expected)
    # Test using series accessor
    expected = pd.Series(result)
    result_series = pd.Series(data_for_encoding()).genomics.encode_additive()
    assert_series_equal(result_series, expected)
    # Test using DataFrame accessor
    df = pd.DataFrame.from_dict(
        {n: data_for_encoding() for n in "ABC"}, orient="columns"
    )
    expected = pd.DataFrame.from_dict({n: result_series for n in "ABC"})
    result_df = df.genomics.encode_additive()
    assert_frame_equal(result_df, expected)


def test_encoding_dominant(data_for_encoding):
    # Test arrays directly
    expected = pd.array([0, 1, 1, 1, None], dtype="UInt8")
    result = data_for_encoding().encode_dominant()
    assert_extension_array_equal(result, expected)
    # Test using series accessor
    expected = pd.Series(result)
    result_series = pd.Series(data_for_encoding()).genomics.encode_dominant()
    assert_series_equal(result_series, expected)
    # Test using DataFrame accessor
    df = pd.DataFrame.from_dict(
        {n: data_for_encoding() for n in "ABC"}, orient="columns"
    )
    expected = pd.DataFrame.from_dict({n: result_series for n in "ABC"})
    result_df = df.genomics.encode_dominant()
    assert_frame_equal(result_df, expected)


def test_encoding_recessive(data_for_encoding):
    # Test arrays directly
    expected = pd.array([0, 0, 1, 1, None], dtype="UInt8")
    result = data_for_encoding().encode_recessive()
    assert_extension_array_equal(result, expected)
    # Test using series accessor
    expected = pd.Series(result)
    result_series = pd.Series(data_for_encoding()).genomics.encode_recessive()
    assert_series_equal(result_series, expected)
    # Test using DataFrame accessor
    df = pd.DataFrame.from_dict(
        {n: data_for_encoding() for n in "ABC"}, orient="columns"
    )
    expected = pd.DataFrame.from_dict({n: result_series for n in "ABC"})
    result_df = df.genomics.encode_recessive()
    assert_frame_equal(result_df, expected)


def test_encoding_codominant(data_for_encoding):
    # Test arrays directly
    expected = pd.Categorical(
        ["Ref", "Het", "Hom", "Hom", None],
        categories=["Ref", "Het", "Hom"],
        ordered=True,
    )
    result = data_for_encoding().encode_codominant()
    assert_extension_array_equal(result, expected)
    # Test using series accessor
    expected = pd.Series(result)
    result_series = pd.Series(data_for_encoding()).genomics.encode_codominant()
    assert_series_equal(result_series, expected)
    # Test using DataFrame accessor
    df = pd.DataFrame.from_dict(
        {n: data_for_encoding() for n in "ABC"}, orient="columns"
    )
    expected = pd.DataFrame.from_dict({n: result_series for n in "ABC"})
    result_df = df.genomics.encode_codominant()
    assert_frame_equal(result_df, expected)


@pytest.mark.parametrize(
    "alpha_value,ref_allele,alt_allele,minor_allele_freq,expected",
    [
        (0.25, "A", "T", 0.45, pd.array([0.0, 0.25, 1.0, None, None], dtype="Float64")),
        (0.25, "T", "A", 0.45, pd.array([1.0, 0.25, 0.0, None, None], dtype="Float64")),
        (1.0, "T", "A", 0.45, pd.array([1.0, 1.0, 0.0, None, None], dtype="Float64")),
        pytest.param(
            0.25,
            "A",
            "G",
            0.45,
            None,
            marks=pytest.mark.xfail(
                raises=ValueError, strict=True, reason="Missing Alt Allele"
            ),
        ),
    ],
)
def test_encoding_weighted(
    data_for_encoding, alpha_value, ref_allele, alt_allele, minor_allele_freq, expected
):
    result = data_for_encoding().encode_weighted(
        alpha_value, ref_allele, alt_allele, minor_allele_freq
    )
    assert_extension_array_equal(expected, result)


@pytest.mark.parametrize(
    "encoding_info,expected",
    [
        # Standard working test
        (
            pd.DataFrame(
                {
                    "Variant ID": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                    "Alpha Value": [0.10, 0.20, 0.30, 0.40, 0.50],
                    "Ref Allele": ["A", "B", "C", "D", "E"],
                    "Alt Allele": ["a", "b", "c", "d", "e"],
                    "Minor Allele Frequency": [0.1, 0.1, 0.1, 0.1, 0.1],
                }
            ),
            pd.DataFrame(
                {
                    "var0": [0.0, 0.1, 1.0, None, None],
                    "var1": [0.0, 0.2, 1.0, None, None],
                    "var2": [0.0, 0.3, 1.0, None, None],
                    "var3": [0.0, 0.4, 1.0, None, None],
                    "var4": [0.0, 0.5, 1.0, None, None],
                },
                dtype="Float64",
            ),
        ),
        # Missing one variant, swapping one variant, wrong allele one variant
        (
            pd.DataFrame(
                {
                    "Variant ID": ["rs1", "rs2", "rs4", "rs5"],
                    "Alpha Value": [0.10, 0.20, 0.40, 0.50],
                    "Ref Allele": ["A", "B", "D", "e"],
                    "Alt Allele": ["a", "b", "X", "E"],
                    "Minor Allele Frequency": [0.1, 0.1, 0.1, 0.1],
                }
            ),
            pd.DataFrame(
                {
                    "var0": [0.0, 0.1, 1.0, None, None],
                    "var1": [0.0, 0.2, 1.0, None, None],
                    "var4": [1.0, 0.5, 0.0, None, None],
                },
                dtype="Float64",
            ),
        ),
    ],
)
def test_encoding_weighted_df(encoding_df, encoding_info, expected):
    result = encoding_df.genomics.encode_weighted(encoding_info)
    assert_frame_equal(expected, result)


def test_generated_encodings_plink(genotypearray_df):
    data = pd.DataFrame(
        {"phenotype": genotypearray_df.index.get_level_values("phenotype")},
        index=genotypearray_df.index,
    )
    result_df = genotypearray_df.genomics.generate_weighted_encodings(
        data, outcome_variable="phenotype"
    )
    result_series = genotypearray_df[
        "18_nullA_18"
    ].genomics.generate_weighted_encodings(data, outcome_variable="phenotype")
    expected = pd.DataFrame(
        {
            "Variant ID": ["nullA_18"],
            "Alpha Value": [0.0],
            "Ref Allele": ["D"],
            "Alt Allele": ["d"],
            "Minor Allele Frequency": [40 / 3000],
        }
    )
    assert_frame_equal(expected, result_df)
    assert_frame_equal(expected, result_series)


@pytest.mark.parametrize(
    "bam,expected_alphas",
    [
        (
            BAMS.from_model(
                SNPEffectEncodings.DOMINANT,
                SNPEffectEncodings.DOMINANT,
                main1=1,
                main2=1,
                interaction=0,
            ),
            [0.884755, 1.292038],
        ),
        (
            BAMS.from_model(
                SNPEffectEncodings.SUPER_ADDITIVE,
                SNPEffectEncodings.SUPER_ADDITIVE,
                main1=1,
                main2=1,
                interaction=0,
            ),
            [0.504485, 0.757154],
        ),
        (
            BAMS.from_model(
                SNPEffectEncodings.ADDITIVE,
                SNPEffectEncodings.ADDITIVE,
                main1=1,
                main2=1,
                interaction=0,
            ),
            [0.21231, 0.41846],
        ),
        (
            BAMS.from_model(
                SNPEffectEncodings.SUB_ADDITIVE,
                SNPEffectEncodings.SUB_ADDITIVE,
                main1=1,
                main2=1,
                interaction=0,
            ),
            [0.091149, 0.163882],
        ),
        (
            BAMS.from_model(
                SNPEffectEncodings.RECESSIVE,
                SNPEffectEncodings.RECESSIVE,
                main1=1,
                main2=1,
                interaction=0,
            ),
            [-0.332334, -0.252329],
        ),
    ],
)
def test_generated_encodings_bams(bam, expected_alphas):
    genotypes = bam.generate_case_control()
    data = genotypes["Outcome"]
    genotypes = genotypes.drop(columns="Outcome")
    result = genotypes.genomics.generate_weighted_encodings(
        data, outcome_variable="Outcome"
    )
    assert np.isclose(result["Alpha Value"], expected_alphas).all()
