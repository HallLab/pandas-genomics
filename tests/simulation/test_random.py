from pandas_genomics.scalars import Variant
from pandas_genomics import sim


def test():
    var = Variant(chromosome="1", position=123456, ref="T", alt=["A"])
    gta = sim.generate_random_gt(var, alt_allele_freq=0.3)
    var2 = Variant(chromosome="1", position=223456, ref="T", alt=["A", "C"])
    gta_2 = sim.generate_random_gt(var2, alt_allele_freq=[0.25, 0.05])
