"""
Author: Erik Kubaczka
E-Mail: erik.kubaczka@tu-darmstadt.de
"""


import numpy as np

from models.promoter_model import PromoterModel

from models.custom_cache import cache_this


class MomentModel:

    def __init__(self,
                 num_promoter_states,
                 num_trainable_parameters,
                 infinitesimal_generator_function,
                 transcription_rates,
                 mRNA_degradation_rate,
                 translation_rate,
                 protein_degradation_rate,
                 normalize_cycle_time=False,
                 external_concentrations_for_normalization=-1):
        self.promoter_model = PromoterModel(num_states=num_promoter_states,
                                            num_trainable_parameters=num_trainable_parameters,
                                            infinitesimal_generator_function=infinitesimal_generator_function,
                                            transcription_rates=transcription_rates,
                                            normalize_cycle_time=normalize_cycle_time,
                                            external_concentrations_for_normalization=external_concentrations_for_normalization)

        self.mu1_rates = transcription_rates
        self.d1 = mRNA_degradation_rate
        self.mu2 = translation_rate
        self.d2 = protein_degradation_rate

    def E_M(self, external_concentrations):
        E_mu1 = self.promoter_model.E_mu1(external_concentrations)
        E_M = E_mu1 / self.d1
        return E_M

    def E_P(self, external_concentrations):
        E_M = self.E_M(external_concentrations)
        E_P = E_M * self.mu2 / self.d2
        return E_P

    @cache_this
    def func(self, external_concentrations, factor):
        N_states = self.promoter_model.N_states
        LAMBDA = self.promoter_model.propensity_matrix(external_concentrations=external_concentrations)

        M = factor * np.eye(N_states) - np.transpose(LAMBDA)
        # M = M.astype(np.longdouble)
        try:
            M_inv = np.linalg.inv(M)
        except np.linalg.LinAlgError:
            print("LinAlgError Occured. Performing Crude Approximation")
            print(M)
            M_inv = np.linalg.inv(M + np.max(M) * np.eye(N_states))
        return M_inv

    @cache_this
    def E_M_Xx(self, external_concentrations):
        M_inv = self.func(external_concentrations=external_concentrations, factor=self.d1)
        E_mu1x_Xx_squared = self.promoter_model.E_mu1x_Xx_squared(external_concentrations=external_concentrations)
        b = E_mu1x_Xx_squared
        E_M_Xx = M_inv @ b

        # E_M_Xx = np.array([0.01,  0.049, 2.238, 1.139])
        return E_M_Xx

    def E_M_mu1(self, external_concentrations):
        vals = self.mu1_rates * self.E_M_Xx(external_concentrations)
        E_M_mu1 = np.sum(vals)
        return E_M_mu1

    def E_P_Xx(self, external_concentrations):
        M_inv = self.func(external_concentrations=external_concentrations, factor=self.d2)

        E_M_Xx = self.E_M_Xx(external_concentrations)

        b = E_M_Xx * 1
        E_P_Xx = self.mu2 * M_inv @ b

        return E_P_Xx

    def E_P_mu1(self, external_concentrations):
        E_P_Xx = self.E_P_Xx(external_concentrations)
        vals = self.mu1_rates * E_P_Xx
        E_P_mu1 = np.sum(vals)

        return E_P_mu1

    def E_M_P(self, external_concentrations):
        E_P_mu1 = self.E_P_mu1(external_concentrations)
        E_M_squared = self.E_M_squared(external_concentrations)

        E_M_P = E_P_mu1 + self.mu2 * E_M_squared
        E_M_P = E_M_P / (self.d1 + self.d2)
        return E_M_P

    @cache_this
    def E_M_squared(self, external_concentrations):
        E_mu1 = self.promoter_model.E_mu1(external_concentrations)
        E_M_mu1 = self.E_M_mu1(external_concentrations)
        E_M = self.E_M(external_concentrations)

        E_M_squared = E_mu1 + 2 * E_M_mu1 + self.d1 * E_M
        E_M_squared = E_M_squared / (2 * self.d1)

        return E_M_squared

    def E_P_squared(self, external_concentrations):
        E_M = self.E_M(external_concentrations)
        E_M_P = self.E_M_P(external_concentrations)
        E_P = self.E_P(external_concentrations)

        E_P_Squared = self.mu2 * E_M + 2 * self.mu2 * E_M_P + self.d2 * E_P
        E_P_Squared = E_P_Squared / (2 * self.d2)

        return E_P_Squared

    def mean_M(self, external_concentrations):
        return self.E_M(external_concentrations)

    def var_M(self, external_concentrations):
        return self.E_M_squared(external_concentrations) - (self.E_M(external_concentrations)) ** 2

    def mean_P(self, external_concentrations):
        return self.E_P(external_concentrations)

    def var_P(self, external_concentrations):
        E_P_squared = self.E_P_squared(external_concentrations)
        squared_E_P = (self.E_P(external_concentrations)) ** 2
        variance = E_P_squared - squared_E_P
        return variance
