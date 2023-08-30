"""
Author: Erik Kubaczka
E-Mail: erik.kubaczka@tu-darmstadt.de
"""

import numpy as np

from models.moment_model import MomentModel


class EnergyAwareMomentModel(MomentModel):
    def __init__(self,
                 gate_parameters,
                 num_trainable_parameters,
                 infinitesimal_generator_function):
        self.gate_parameters = gate_parameters
        N_promoter_states = gate_parameters["GENE"]["N_PROMOTER_STATES"]
        transcription_rates = gate_parameters["mRNA"]["TRANSCRIPTION_RATES"]
        mRNA_degradation_rate = gate_parameters["mRNA"]["mRNA_DEGRADATION_RATE"]
        translation_rate = gate_parameters["PROTEIN"]["TRANSLATION_RATE"]
        protein_degradation_rate = gate_parameters["PROTEIN"]["PROTEIN_DEGRADATION_RATE"]

        super().__init__(N_promoter_states,
                         num_trainable_parameters=num_trainable_parameters,
                         infinitesimal_generator_function=infinitesimal_generator_function,
                         transcription_rates=transcription_rates,
                         mRNA_degradation_rate=mRNA_degradation_rate,
                         translation_rate=translation_rate,
                         protein_degradation_rate=protein_degradation_rate,
                         normalize_cycle_time=False,
                         external_concentrations_for_normalization=-1)

        self.gate_parameters = gate_parameters

        # Energy data expects data of the form
        # {"mRNA": {
        #     "TRANSCRIPTION_RATES": [
        #         0.0,
        #         0.0,
        #         100000.0,
        #         100000.0
        #     ],
        #     "mRNA_DEGRADATION_RATE": 1,
        #     "e_m": 16,
        #     "e_m_const": 0
        # },
        #     "PROTEIN": {
        #         "TRANSLATION_RATE": 1,
        #         "PROTEIN_DEGRADATION_RATE": 1,
        #         "e_p": 42,
        #         "e_p_const": 52
        #     }
        # }

        pass

    def energy_rate(self, external_concentrations):
        mean_P, var_P, mean_M, var_M, energy_rate = self.get_distributions_and_energy_rate(self, external_concentrations)
        return energy_rate


    def get_distributions_and_energy_rate(self, external_concentrations):
        def get_energy_per_molecule(molecule_data):
            if "l" not in molecule_data:
                print(r"The molecule length is not povided. Setting $l=0$.")
                l = 0
            else:
                l = molecule_data["l"]

            if "e" not in molecule_data:
                print(r"The per element energy consumption is not provided. Setting $e=0$.")
                e = 0
            else:
                e = molecule_data["e"]

            if "e_const" not in molecule_data:
                print(r"The constant per molecule energy consumption is not provided. Setting $e_{const}=0$.")
                e_const = 0
            else:
                e_const = molecule_data["e_const"]

            return e * l + e_const

        if self.gate_parameters is None:
            return np.nan

        mRNA_mean = self.mean_M(external_concentrations=external_concentrations)
        mRNA_variance = self.var_M(external_concentrations=external_concentrations)

        protein_mean = self.mean_P(external_concentrations=external_concentrations)
        protein_variance = self.var_P(external_concentrations=external_concentrations)

        epsilon_p = self.promoter_model.entropy_production_rate(external_concentrations=external_concentrations)
        e_mRNA = get_energy_per_molecule(self.gate_parameters["mRNA"])
        e_protein = get_energy_per_molecule(self.gate_parameters["PROTEIN"])

        mRNA_degradation_rate = self.d1
        protein_degradation_rate = self.d2
        e_tx = e_mRNA * mRNA_degradation_rate * mRNA_mean
        e_tl = e_protein * protein_degradation_rate * protein_mean

        cur_energy_rate = epsilon_p + e_tx + e_tl

        return protein_mean, protein_variance, mRNA_mean, mRNA_variance, cur_energy_rate