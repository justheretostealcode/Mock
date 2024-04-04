"""
Author: Erik Kubaczka
E-Mail: erik.kubaczka@tu-darmstadt.de
"""
import numpy as np
from scipy import linalg


class RNAMomentModel:
    def __init__(self, transcription_rate, degradation_rate, energy_per_nucleotide, energy_per_rna, length):
        """
        :param transcription_rate: The transcription rate in RNA per second
        :param degradation_rate: The degradation rate in RNA per second
        :param energy_per_nucleotide: The energy required to extend the nucleotide chain by a single nucleotide or the length dependent part of RNA degradation
        :param energy_per_rna: The length independent energy required for a single RNA during translation or degradation
        :param length: The length of the RNA chain in nucleotides
        """
        self.transcription_rate = transcription_rate
        self.degradation_rate = degradation_rate
        self.energy_per_nucleotide = energy_per_nucleotide
        self.energy_per_rna = energy_per_rna
        self.length = length
        self.energy_per_rna = (self.energy_per_nucleotide * self.length + self.energy_per_rna)

    def __call__(self, promoter_output: dict, sim_settings: dict, *args, **kwargs) -> dict:
        mode = "det"
        if "mode" in sim_settings:
            mode = sim_settings["mode"]

        rna_output = dict(promoter_output)

        propensity_matrix = promoter_output["propensity_matrix"]
        avg_promoter_activity = promoter_output["average_promoter_activity"]
        avg_promoter_activity_per_state = promoter_output["average_promoter_activity_per_state"]
        promoter_activity_per_state = promoter_output["promoter_activity_per_state"]

        n_samples = avg_promoter_activity.shape[0]

        l = self.transcription_rate / self.degradation_rate

        """
        Deriving the mean RNA
        """
        self.mean = l * avg_promoter_activity

        """
        Deriving the variance of the RNA
        """
        self.var = np.nan * np.empty(shape=self.mean.shape)
        M_inv = None
        if mode != "det":
            M = self.degradation_rate * np.eye(*propensity_matrix.shape[1:]) - propensity_matrix.transpose((0, 2, 1))
            # M_inv = np.linalg.inv(M)
            # rna_promoter_state_correlation_legacy = M_inv @ np.expand_dims(avg_promoter_activity_per_state, axis=-1)
            # rna_promoter_state_correlation_legacy = self.transcription_rate * rna_promoter_state_correlation_legacy.squeeze(axis=-1)
            rna_promoter_state_correlation = np.linalg.solve(M, avg_promoter_activity_per_state)
            rna_promoter_state_correlation = self.transcription_rate * rna_promoter_state_correlation
            # Equation is checked
            E_rna_squared = np.expand_dims(promoter_activity_per_state, axis=-2) @ np.expand_dims(
                rna_promoter_state_correlation, axis=-1)
            E_rna_squared = self.mean + l * E_rna_squared[:, 0, 0]  # Performs squeezing of axis 1 and 2

            self.var = E_rna_squared - self.mean ** 2

        """
        Energy Considerations
        """
        rna_production_and_degradation_rate = self.degradation_rate * self.mean

        self.mean_energy_dissipation_rate = rna_production_and_degradation_rate * self.energy_per_rna

        rna_output["rna_mean"] = self.mean
        rna_output["rna_var"] = self.var
        rna_output["M"] = M
        rna_output["transcription_rate"] = self.transcription_rate
        rna_output["rna_degradation_rate"] = self.degradation_rate
        rna_output["rna_mean_energy_dissipation_rate"] = self.mean_energy_dissipation_rate
        rna_output["energy_tx"] = self.mean_energy_dissipation_rate

        return rna_output


class ProteinMomentModel:
    def __init__(self, translation_rate, degradation_rate, energy_per_amino_acid, energy_per_protein, length):
        """

        :param translation_rate: The translation rate in proteins per second
        :param degradation_rate: The degradation rate in proteins per second
        :param energy_per_amino_acid: The energy required to extend the polypeptide chain by a single amino acid or the length dependent part of protein degradation
        :param energy_per_protein: The length independent energy required for a single protein during translation or degradation
        :param length: The length of the polypeptide chain in amino acids
        """
        self.translation_rate = translation_rate
        self.degradation_rate = degradation_rate
        self.energy_per_amino_acid = energy_per_amino_acid
        self.energy_per_protein = energy_per_protein
        self.length = length

        self.energy_per_protein = self.energy_per_amino_acid * self.length + self.energy_per_protein

    def __call__(self, rna_output: dict, sim_settings: dict, *args, **kwargs) -> dict:
        mode = "det"
        if "mode" in sim_settings:
            mode = sim_settings["mode"]

        propensity_matrix = rna_output["propensity_matrix"]
        avg_promoter_activity_per_state = rna_output["average_promoter_activity_per_state"]
        promoter_activity_per_state = rna_output["promoter_activity_per_state"]
        transcription_rate = rna_output["transcription_rate"]
        rna_degradation_rate = rna_output["rna_degradation_rate"]
        M = rna_output["M"]
        rna_mean = rna_output["rna_mean"]
        rna_var = rna_output["rna_var"]

        l = self.translation_rate / self.degradation_rate

        """
        Derivation of protein's mean
        """
        self.mean = l * rna_mean

        """
        Derivation of protein's variance
        """

        self.var = np.nan * np.ones(shape=self.mean.shape)
        self.protein_rna_covariance = np.nan * np.ones(shape=self.mean.shape)
        # raise Exception("Reimplement on the basis of E[n_m * n_p]")

        if mode != "det" and M is not None:
            # M_inv = np.linalg.inv(M)
            M_2 = self.degradation_rate * np.eye(*propensity_matrix.shape[1:]) - propensity_matrix.transpose((0, 2, 1))
            # M_inv_2 = np.linalg.inv(M_2)

            denominator = rna_degradation_rate + self.degradation_rate

            # Transforms it into a matrix representing a row vector
            # promoter_activity_per_state_legacy = np.expand_dims(promoter_activity_per_state, axis=-2)
            # Transforms it into a matrix representing a column vector
            # avg_promoter_activity_per_state_legacy = np.expand_dims(avg_promoter_activity_per_state, axis=-1)

            # One could save a matrix multiplication by reusing a result from the RNAMomentModel
            # protein_promoter_correlation = self.translation_rate * transcription_rate * M_inv_2 @ M_inv \ @ avg_promoter_activity_per_state
            # protein_promoter_correlation_legacy = self.translation_rate * transcription_rate * M_inv_2 @ M_inv @ avg_promoter_activity_per_state_legacy
            protein_promoter_correlation = np.linalg.solve(M @ M_2, self.translation_rate * transcription_rate * avg_promoter_activity_per_state)
            E_m_squared = rna_var + rna_mean ** 2

            E_m_p_legacy = self.translation_rate * E_m_squared
            E_m_p = self.translation_rate * E_m_squared
            # [:, 0, 0] Squeezes axis -1 and -2
            # E_m_p_legacy += transcription_rate * (promoter_activity_per_state_legacy @ protein_promoter_correlation_legacy)[:, 0, 0]
            E_m_p += transcription_rate * (np.expand_dims(promoter_activity_per_state, axis=-2) @ np.expand_dims(protein_promoter_correlation, axis=-1))[:, 0, 0]
            # E_m_p_legacy = E_m_p_legacy / denominator
            E_m_p = E_m_p / denominator

            E_p_squared = l * E_m_p + self.mean

            self.var = E_p_squared - self.mean ** 2
            self.protein_rna_covariance = E_m_p - rna_mean * self.mean
            # self.var = self.mean
            # self.var -= self.mean ** 2
            # self.var += self.translation_rate / denominator * self.mean
            # self.var += self.translation_rate * l * transcription_rate ** 2 / denominator * modified_protein_promoter_correlation
            #
            # self.protein_rna_covariance = (self.translation_rate * (rna_var + rna_mean ** 2)
            #                                + transcription_rate * promoter_activity_per_state @ protein_promoter_correlation) * denominator
            # self.protein_rna_covariance -= self.mean * rna_mean

            self.var = self.var
            self.protein_rna_covariance = self.protein_rna_covariance

        """
        Energy Considerations
        """
        protein_production_and_degradation_rate = self.degradation_rate * self.mean

        self.mean_energy_dissipation_rate = protein_production_and_degradation_rate * self.energy_per_protein

        protein_output = dict(rna_output)
        protein_output["protein_mean"] = self.mean
        protein_output["protein_var"] = self.var
        protein_output["protein_rna_covariance"] = self.protein_rna_covariance
        protein_output["protein_mean_energy_dissipation_rate"] = self.mean_energy_dissipation_rate
        protein_output["energy_tl"] = self.mean_energy_dissipation_rate

        return protein_output


class CombinedMomentModel:
    def __init__(self, transcription_rate, rna_degradation_rate, energy_per_nucleotide, energy_per_rna, rna_length,
                 translation_rate, protein_degradation_rate, energy_per_amino_acid, energy_per_protein, protein_length):
        self.rna_model = RNAMomentModel(transcription_rate, rna_degradation_rate, energy_per_nucleotide, energy_per_rna,
                                        rna_length)
        self.protein_model = ProteinMomentModel(translation_rate, protein_degradation_rate, energy_per_amino_acid,
                                                energy_per_protein, protein_length)

        l_m = self.rna_model.transcription_rate / self.rna_model.degradation_rate
        l_p = self.protein_model.translation_rate / self.protein_model.degradation_rate
        self.scaling_factor = l_m * l_p  # Factor for converting RPU to output
        pass

    def __call__(self, promoter_output: dict, sim_settings: dict, *args, **kwargs) -> dict:
        rna_output = self.rna_model(promoter_output, sim_settings=sim_settings)
        protein_output = self.protein_model(rna_output, sim_settings=sim_settings)
        # raise Exception("Check correctness of outcome!")
        return protein_output
