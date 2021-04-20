import numpy as np

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
    simulated_df = bas.generate_case_control()
    # Count unique rows
    counts = (
        simulated_df.groupby(["Outcome", "SNP1", "SNP2"])
        .size()
        .reset_index(name="Count")
    )
    print()
