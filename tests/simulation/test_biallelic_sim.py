from pandas_genomics import BAMS, SNPEffectEncodings, PenetranceTables


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
    print(bas)
    test_sim = BAMS.from_model(
        SNPEffectEncodings.RECESSIVE,
        SNPEffectEncodings.RECESSIVE,
        main1=1,
        main2=1,
        interaction=1,
    )
    simulated_df_cc = test_sim.generate_case_control(snr=0.99)
    simulated_df_quant = test_sim.generate_quantitative()
