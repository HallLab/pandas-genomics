"""
Simulation
----------
The `sim` module provides classes for generating simulated genotypes

    .. autosummary::
        :toctree: sim

        BAMS
        SNPEffectEncodings
        PenetranceTables
        generate_random_gt

"""


from .biallelic_model_simulator import BAMS, SNPEffectEncodings, PenetranceTables
from .random_gt import generate_random_gt
