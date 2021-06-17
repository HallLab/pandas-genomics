from pandas._testing import assert_frame_equal

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


def test():
    assert BAMS() == BAMS.from_model()
    assert BAMS(PenetranceTables.HET_HET) == BAMS.from_model(
        eff1=(0, 1, 0), eff2=(0, 1, 0), main1=0, main2=0, interaction=1
    )
    bas = BAMS.from_model(
        SNPEffectEncodings.DOMINANT,
        SNPEffectEncodings.DOMINANT,
        main1=0,
        main2=0,
        interaction=1,
    )
    test_sim = BAMS.from_model(
        SNPEffectEncodings.RECESSIVE,
        SNPEffectEncodings.RECESSIVE,
        main1=1,
        main2=1,
        interaction=1,
    )

    # Test simulating data using random seeds
    simulated_df_cc = test_sim.generate_case_control(snr=0.01)
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
