import numpy as np

from pandas_genomics import BialleleicSimulation, SNPEffectEncodings


def test():
    bas = BialleleicSimulation()
    print(bas)
    print()
    bas2 = BialleleicSimulation.from_model(SNPEffectEncodings.DOMINANT,
                                           SNPEffectEncodings.DOMINANT,
                                           main1=5, main2=5, interaction=0)
    print(bas2)
    simulated_df = bas2.generate_case_control()
    # Count unique rows
    counts = simulated_df.groupby(['Outcome', 'SNP1', 'SNP2']).size().reset_index(name='Count')
    print()

