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
    bas2.generate_case_control()

