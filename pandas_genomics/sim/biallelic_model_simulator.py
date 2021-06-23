from enum import Enum
from typing import Optional, Tuple, Union

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm

from pandas_genomics.arrays import GenotypeArray, GenotypeDtype
from pandas_genomics.scalars import Variant


class SNPEffectEncodings(Enum):
    """Enum: Normalized SNP Effects encoded as 3-length tuples"""

    DOMINANT = (0, 1, 1)
    SUPER_ADDITIVE = (0, 0.75, 1)
    ADDITIVE = (0, 0.5, 1)
    SUB_ADDITIVE = (0, 0.25, 1)
    RECESSIVE = (0, 0, 1)
    HET = (0, 1, 0)


class PenetranceTables(Enum):
    """Enum: Penetrance Tables for Simple Models"""

    HR_HR = [1, 0, 0, 0, 0, 0, 0, 0, 0]  # Homozygous Referent X Homozygous Referent
    HR_HET = [0, 1, 0, 0, 0, 0, 0, 0, 0]  # Homozygous Referent X Heterozygous
    HR_HA = [0, 0, 1, 0, 0, 0, 0, 0, 0]  # Homozygous Referent X Homozygous Alternate
    HET_HET = [0, 0, 0, 0, 1, 0, 0, 0, 0]  # Heterozygous X Heterozygous
    HET_HA = [0, 0, 0, 0, 0, 1, 0, 0, 0]  # Heterozygous X Homozygous Alternate
    XOR = [1, 0, 1, 0, 1, 0, 1, 0, 1]  # XOR Model
    HYP = [0, 0.5, 1, 0.5, 0.5, 0.5, 1, 0.5, 0]  # Hyperbolic Model
    RHYP = [1, 0.5, 0, 0.5, 0.5, 0.5, 0, 0.5, 1]  # Hyperbolic Model
    NULL = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]  # Null Model (Always 50/50)


class BAMS:
    """
    Biallelic Model Simulator.  Used to simulate two SNPs with phenotype data based on a penetrance table.

    It can be initialized using the PenetranceTables enum or using `from_model` with values from the SNPEffectEncodings enum.
    """

    def __init__(
        self,
        pen_table: Union[np.array, PenetranceTables] = np.array(
            [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 1.0, 2.0]]
        ),
        snp1: Optional[Variant] = None,
        snp2: Optional[Variant] = None,
        random_seed: int = 1855,
    ):
        pen_table, snp1, snp2 = self._validate_params(pen_table, snp1, snp2)
        self.pen_table = pen_table
        self.snp1 = snp1
        self.snp2 = snp2
        self._random_seed = random_seed

    def __str__(self):
        pen_table_df = pd.DataFrame(self.pen_table)
        pen_table_df.columns = self._get_genotype_strs(self.snp1)
        pen_table_df.index = self._get_genotype_strs(self.snp2)
        return (
            f"SNP1 = {str(self.snp1)}\n"
            f"SNP2 = {str(self.snp2)}\n"
            f"Penetrance Table:\n"
            f"----------------\n"
            f"{pen_table_df}\n"
            f"----------------\n"
            f"Random Seed = {self._random_seed}"
        )

    def __eq__(self, other):
        if type(other) is not type(self):
            raise NotImplemented
        else:
            return (
                (self.pen_table == other.pen_table).all()
                & (self.snp1 == other.snp1)
                & (self.snp2 == other.snp2)
                & (self._random_seed == other._random_seed)
            )

    @property
    def random_seed(self) -> int:
        return self._random_seed

    def set_random_seed(self, new_seed: int):
        self._random_seed = new_seed

    @staticmethod
    def _validate_params(pen_table, snp1, snp2):
        # Process Enum
        if type(pen_table) is PenetranceTables:
            pen_table = np.array(pen_table.value).reshape((3, 3))
        elif pen_table.shape != (3, 3):
            raise ValueError(f"Incorrect shape for pen_table, must be 3x3")

        # SNPs
        if snp1 is None:
            snp1 = Variant(id="rs1", ref="A", alt=["a"])
        if snp2 is None:
            snp2 = Variant(id="rs2", ref="B", alt=["b"])

        if len(snp1.alt) != 1:
            raise ValueError(f"SNP1 is not Bialleleic: {snp1}")
        if len(snp2.alt) != 1:
            raise ValueError(f"SNP2 is not Bialleleic: {snp2}")

        return pen_table, snp1, snp2

    @classmethod
    def from_model(
        cls,
        eff1: Union[
            Tuple[float, float, float], SNPEffectEncodings
        ] = SNPEffectEncodings.RECESSIVE,
        eff2: Union[
            Tuple[float, float, float], SNPEffectEncodings
        ] = SNPEffectEncodings.RECESSIVE,
        baseline: float = 0.0,
        main1: float = 1.0,
        main2: float = 1.0,
        interaction: float = 0.0,
        snp1: Optional[Variant] = None,
        snp2: Optional[Variant] = None,
        random_seed: int = 1855,
    ):
        """
        Create a BiallelicSimulation with a Penetrance Table based on a fully specified model
        y = β0 + β1(eff1) + β2(eff2) + β3(eff1*eff2)

        Parameters
        ----------
        eff1: tuple of 3 floats
            Normalized effect of SNP1
            May be passed a value from the `Effects` enum
        eff2: tuple of 3 floats
            Normalized effect of SNP2
            May be passed a value from the `Effects` enum
        baseline: float, default 0.0
            β0 in the formula
        main1: float, default 1.0
            Main effect of SNP1, β1 in the formula
        main2: float, default 1.0
            Main effect of SNP2, β2 in the formula
        interaction: float, default 0.0
            Interaction effect, β3 in the formula
        snp1: Variant, default is None
            First SNP, one will be created by default if not specified
        snp2: Variant, default is None
            Second SNP, one will be created by default if not specified
        random_seed: int, default is 1855
            Random seed used during simulation


        Returns
        -------
        BAMS

        """
        # TODO: Add more validation
        if type(eff1) is SNPEffectEncodings:
            eff1 = eff1.value
        if type(eff2) is SNPEffectEncodings:
            eff2 = eff2.value

        # Shape effects and scale if needed
        eff1 = np.array([eff1])  # SNP1 = columns
        eff2 = np.array([eff2]).transpose()  # SNP2 = rows
        if eff1.min() != 0 or eff1.max() != 1:
            print("Scaling eff1")
            eff1 = (eff1 - eff1.min()) / (eff1.max() - eff1.min())
        if eff2.min() != 0 or eff2.max() != 1:
            print("Scaling eff2")
            eff2 = (eff2 - eff2.min()) / (eff2.max() - eff2.min())

        pen_table = (
            baseline
            + main1 * np.repeat(eff1, 3, axis=0)
            + main2 * np.repeat(eff2, 3, axis=1)
            + interaction * np.outer(eff2, eff1)
        )

        return cls(pen_table, snp1, snp2, random_seed)

    def generate_case_control(
        self,
        n_cases: int = 1000,
        n_controls: int = 1000,
        maf1: float = 0.30,
        maf2: float = 0.30,
        snr: Optional[float] = None,
    ):
        """
        Simulate genotypes with the specified number of 'case' and 'control' phenotypes

        Parameters
        ----------
        n_cases: int, default 1000
        n_controls: int, default 1000
        maf1: float, default 0.30
            Minor Allele Frequency to use for SNP1
        maf2: float, default 0.30
            Minor Allele Frequency to use for SNP2
        snr: float, default 1.0
            Signal-to-noise ratio

        Returns
        -------
        pd.Dataframe
            Dataframe with 3 columns: Outcome (categorical), SNP1 (GenotypeArray), and SNP2 (GenotypeArray)

        """
        np.random.seed(self._random_seed)
        # Validate params
        if n_cases < 1 or n_controls < 0:
            raise ValueError(
                "Simulation must include at least one case and at least one control"
            )

        # Scale the penetrance table to be between 0 and 1
        pen_table = self.pen_table
        if (pen_table.min() != 0) or (pen_table.max() != 1):
            pen_table_min = pen_table.min()
            pen_table_range = pen_table.max() - pen_table_min
            if pen_table_range > 0:
                pen_table = (pen_table - pen_table_min) / pen_table_range
                # Otherwise the penetrance table is flat, i.e. a null model

        # Create table of Prob(GT) based on MAF, assuming HWE
        prob_snp1 = np.array([(1 - maf1) ** 2, 2 * maf1 * (1 - maf1), (maf1) ** 2])
        prob_snp2 = np.array(
            [(1 - maf2) ** 2, 2 * maf2 * (1 - maf2), (maf2) ** 2]
        ).transpose()
        prob_gt = np.outer(prob_snp2, prob_snp1)

        if snr is not None:
            sigma = self._calculate_sigma(pen_table, prob_gt)
            # Scale the penetrance by the amount of unexplained variance
            # Odds Ratio = exp(SNP/Sigma)
            min_p = 1 / (1 + np.exp(1 / sigma * snr))
            p_diff = 1 - 2 * min_p
            # Adjust the penetrance
            pen_table = min_p + pen_table * p_diff

        #                     P(Case|GT) * P(GT)
        # Bayes: P(GT|Case) = ------------------
        #                           P(Case)

        # Prob(Case|GT) = pen_table
        # Prob(Case) = sum(Prob(Case|GTi) * Prob(GTi) for each GT i)
        prob_case = (pen_table * prob_gt).sum()
        # Prob(GT|Case)
        prob_gt_given_case = (pen_table * prob_gt) / prob_case

        # Prob(Control|GT) = 1-pen_table
        # Prob(Control) = sum(Prob(Control|GTi) * Prob(GTi) for each GT i)
        prob_control = ((1 - pen_table) * prob_gt).sum()
        # Prob(GT|Control)
        prob_gt_given_control = ((1 - pen_table) * prob_gt) / prob_control

        # Generate genotypes based on the simulated cases and controls
        # Pick int index into the table (0 through 8) counted left to right then top to bottom (due to flatten())
        # TODO: This isn't always consistent in multiple runs
        case_gt_table_idxs = np.random.choice(
            range(9), size=n_cases, p=prob_gt_given_case.flatten()
        )
        control_gt_table_idxs = np.random.choice(
            range(9), size=n_controls, p=prob_gt_given_control.flatten()
        )

        # Create GenotypeArrays
        snp1_case_array = pd.Series(self._get_snp1_gt_array(case_gt_table_idxs))
        snp2_case_array = pd.Series(self._get_snp2_gt_array(case_gt_table_idxs))
        snp1_control_array = pd.Series(self._get_snp1_gt_array(control_gt_table_idxs))
        snp2_control_array = pd.Series(self._get_snp2_gt_array(control_gt_table_idxs))

        # Merge data together
        snp1 = pd.concat([snp1_case_array, snp1_control_array]).reset_index(drop=True)
        snp2 = pd.concat([snp2_case_array, snp2_control_array]).reset_index(drop=True)

        # Generate outcome
        outcome = pd.Series(["Case"] * n_cases + ["Control"] * n_controls).astype(
            "category"
        )
        result = pd.concat([outcome, snp1, snp2], axis=1)
        result.columns = ["Outcome", "SNP1", "SNP2"]

        # Scramble outcome so cases and controls are mixed
        result = result.sample(frac=1, random_state=self._random_seed).reset_index(
            drop=True
        )

        return result

    def generate_quantitative(
        self,
        n_samples: int = 1000,
        maf1: float = 0.30,
        maf2: float = 0.30,
        snr: Optional[float] = None,
    ):
        """
        Simulate genotypes with the specified number of 'case' and 'control' phenotypes

        Parameters
        ----------
        n_samples: int, default 1000
        maf1: float, default 0.30
            Minor Allele Frequency to use for SNP1
        maf2: float, default 0.30
            Minor Allele Frequency to use for SNP2
        snr: float, default 1.0
            Signal-to-noise ratio

        Returns
        -------
        pd.Dataframe
            Dataframe with 3 columns: Outcome, SNP1 (GenotypeArray), and SNP2 (GenotypeArray)

        """
        np.random.seed(self._random_seed)
        pen_table = self.pen_table

        # Create table of Prob(GT) based on MAF, assuming HWE
        prob_snp1 = np.array([(1 - maf1) ** 2, 2 * maf1 * (1 - maf1), (maf1) ** 2])
        prob_snp2 = np.array(
            [(1 - maf2) ** 2, 2 * maf2 * (1 - maf2), (maf2) ** 2]
        ).transpose()
        prob_gt = np.outer(prob_snp2, prob_snp1)

        if snr is not None:
            sigma = self._calculate_sigma(pen_table, prob_gt)
            # Scale the penetrance by the amount of unexplained variance
            pen_table = pen_table / (sigma * snr)

        # Generate genotypes
        gt_table_idxs = np.random.choice(range(9), size=n_samples, p=prob_gt.flatten())
        snp1_array = pd.Series(self._get_snp1_gt_array(gt_table_idxs))
        snp2_array = pd.Series(self._get_snp2_gt_array(gt_table_idxs))

        # Calculate outcome
        outcome = pd.Series(np.take(pen_table.flatten(), gt_table_idxs))

        result = pd.concat([outcome, snp1_array, snp2_array], axis=1)
        result.columns = ["Outcome", "SNP1", "SNP2"]

        return result

    @staticmethod
    def _calculate_sigma(pen_table, prob_gt):
        # Set up a dataframe with the penetrance table outcome for each SNP combination (Encoded as codominant)
        cats = ["Ref", "Het", "Hom"]
        d = pd.DataFrame(
            {
                "SNP1": pd.Categorical(cats * 3, categories=cats, ordered=True),
                "SNP2": pd.Categorical(
                    [i for c in cats for i in [c] * 3],
                    categories=cats,
                    ordered=True,
                ),
            },
            dtype="category",
        )
        d["y"] = pen_table.flatten()
        # Get the weight vec as variance based on genotype probability
        wt_vec = (prob_gt * (1 - prob_gt)).flatten()
        # Regress the possible combinations of genotypes against the penetrance table, using the genotype frequency variance as weights
        y, X = patsy.dmatrices("y ~ 1 + SNP1 + SNP2", data=d)
        mod = sm.WLS(y, X, weights=wt_vec)
        res = mod.fit()
        sigma = np.sqrt(res.scale)  # sigma is sqrt of the dispersion (aka scale)
        return sigma

    @staticmethod
    def _get_genotype_strs(variant):
        """Return a list of homozygous-ref, het, and homozygous-alt"""
        return [
            f"{variant.ref}{variant.ref}",
            f"{variant.ref}{variant.alt[0]}",
            f"{variant.alt[0]}{variant.alt[0]}",
        ]

    def _get_snp1_gt_array(self, gt_table_idxs):
        """Assemble a GenotypeArray for SNP1 directly from genotype table indices"""
        dtype = GenotypeDtype(self.snp1)
        gt_table_data = (
            ((0, 0), np.nan),
            ((0, 1), np.nan),
            ((1, 1), np.nan),
            ((0, 0), np.nan),
            ((0, 1), np.nan),
            ((1, 1), np.nan),
            ((0, 0), np.nan),
            ((0, 1), np.nan),
            ((1, 1), np.nan),
        )
        data = np.array(
            [gt_table_data[i] for i in gt_table_idxs], dtype=dtype._record_type
        )
        return GenotypeArray(values=data, dtype=dtype)

    def _get_snp2_gt_array(self, gt_table_idxs):
        """Assemble a GenotypeArray for SNP2 directly from genotype table indices"""
        dtype = GenotypeDtype(self.snp2)
        gt_table_data = (
            ((0, 0), np.nan),
            ((0, 0), np.nan),
            ((0, 0), np.nan),
            ((0, 1), np.nan),
            ((0, 1), np.nan),
            ((0, 1), np.nan),
            ((1, 1), np.nan),
            ((1, 1), np.nan),
            ((1, 1), np.nan),
        )
        data = np.array(
            [gt_table_data[i] for i in gt_table_idxs], dtype=dtype._record_type
        )
        return GenotypeArray(values=data, dtype=dtype)
