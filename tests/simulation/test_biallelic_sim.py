import numpy as np
import pytest
from pandas._testing import assert_frame_equal

from pandas_genomics import sim
from pandas_genomics.sim import BAMS, SNPEffectEncodings, PenetranceTables


def assert_frame_not_equal(*args, **kwargs):
    try:
        assert_frame_equal(*args, **kwargs)
    except AssertionError:
        # frames are not equal
        pass
    else:
        # frames are equal
        raise AssertionError


@pytest.mark.parametrize(
    "pen_table,baseline,diff,expected",
    [
        (sim.PenetranceTables.NULL, 0.1, 0.8, [0.5] * 9),
        (sim.PenetranceTables.NULL, 0.0, 0.5, [0.25] * 9),
        (
            sim.PenetranceTables.HET_HET,
            0.25,
            0.5,
            [0.25, 0.25, 0.25, 0.25, 0.75, 0.25, 0.25, 0.25, 0.25],
        ),
        (
            np.array(sim.PenetranceTables.HET_HET.value).reshape((3, 3)) * 10,
            0.25,
            0.5,
            [0.25, 0.25, 0.25, 0.25, 0.75, 0.25, 0.25, 0.25, 0.25],
        ),
        (
            sim.PenetranceTables.HET_HA,
            0.1,
            0.9,
            [0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1],
        ),
        pytest.param(
            sim.PenetranceTables.HET_HET,
            -1,
            0,
            [0.0] * 9,
            marks=pytest.mark.xfail(raises=ValueError, strict=True),
        ),
        pytest.param(
            sim.PenetranceTables.HET_HET,
            0.1,
            0.91,
            [0.0] * 9,
            marks=pytest.mark.xfail(raises=ValueError, strict=True),
        ),
        pytest.param(
            [0, 0, 0, 0, 0, 0, 0, 0, 0],
            0.1,
            0.91,
            [0.0] * 9,
            marks=pytest.mark.xfail(raises=ValueError, strict=True),
        ),
    ],
)
def test_pen_table_direct(pen_table, baseline, diff, expected):
    """Test calculation of final penetrance table when a penetrance table is specified"""
    model = BAMS(pen_table=pen_table, penetrance_base=baseline, penetrance_diff=diff)
    np.isclose(model.pen_table, np.reshape(np.array(expected), newshape=(3, 3))).all()


@pytest.mark.parametrize(
    "eff1,eff2,baseline,diff,main1,main2,interaction,expected",
    [
        (
            sim.SNPEffectEncodings.DOMINANT,
            sim.SNPEffectEncodings.DOMINANT,
            0.1,
            0.8,
            1,
            1,
            0,
            [0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.5, 0.9, 0.9],
        ),
        (
            sim.SNPEffectEncodings.DOMINANT,
            sim.SNPEffectEncodings.DOMINANT,
            0.1,
            0.8,
            0,
            0,
            1,
            [0.1, 0.1, 0.1, 0.1, 0.9, 0.9, 0.1, 0.9, 0.9],
        ),
        (
            sim.SNPEffectEncodings.DOMINANT,
            sim.SNPEffectEncodings.DOMINANT,
            0.1,
            0.8,
            1,
            1,
            10,
            [0.1, 1 / 6, 1 / 6, 1 / 6, 0.9, 0.9, 1 / 6, 0.9, 0.9],
        ),
    ],
)
def test_pen_table_model(
    eff1, eff2, baseline, diff, main1, main2, interaction, expected
):
    """Test calculation of final penetrance table from a model"""
    model = BAMS.from_model(
        eff1=eff1,
        eff2=eff2,
        penetrance_base=baseline,
        penetrance_diff=diff,
        main1=main1,
        main2=main2,
        interaction=interaction,
    )
    test = model.generate_quantitative(100, snr=0.1)
    assert np.isclose(
        model.pen_table, np.reshape(np.array(expected), newshape=(3, 3))
    ).all()


def test():
    assert BAMS() == BAMS.from_model()
    assert BAMS(PenetranceTables.HET_HET) == BAMS.from_model(
        eff1=(0, 1, 0), eff2=(0, 1, 0), main1=0, main2=0, interaction=1
    )
    test_sim = BAMS.from_model(
        SNPEffectEncodings.RECESSIVE,
        SNPEffectEncodings.RECESSIVE,
        main1=1,
        main2=1,
        interaction=1,
    )

    # Test simulating data using random seeds
    simulated_df_cc = test_sim.generate_case_control(snr=0.1)
    simulated_df_cc_2 = test_sim.generate_case_control(snr=0.1)
    assert_frame_equal(simulated_df_cc, simulated_df_cc_2)
    test_sim.set_random_seed(123)
    simulated_df_cc_3 = test_sim.generate_case_control(snr=0.1)
    assert_frame_not_equal(simulated_df_cc_2, simulated_df_cc_3)

    # Test quantitative sim
    # TODO: Quantitative model is unfinished
    simulated_df_quant = test_sim.generate_quantitative()


def test_null():
    bas = BAMS(PenetranceTables.NULL)
    simulated = bas.generate_case_control(10000, 1000, 0.1, 0.1)
    # maf should be similar to the specified one despite a large fraction of cases
    # specifically assert it is within 5%
    assert abs(0.1 - simulated["SNP1"].genomics.maf) / 0.1 < 0.05
