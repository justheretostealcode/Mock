"""
Author: Erik Kubaczka
E-Mail: erik.kubaczka@tu-darmstadt.de
"""
import ast
import itertools

import numpy as np

from models.custom_cache import clear_cache, cache_this


class SteadyStateCTMC:

    def __init__(self,
                 N_states,
                 infinitesimal_generator_function):
        self.N_states = N_states
        self.update_infinitesimal_generator_function(infinitesimal_generator_function)
        # self.infinitesimal_generator_function = infinitesimal_generator_function

        self._cur_in_val_dict = None
        self.propensity_matrix = None
        self.distribution = None
        self.entropy_production_rate = None
        pass

    def update_infinitesimal_generator_function(self, new_function: dict):
        self.infinitesimal_generator_function = new_function

        self.infinitesimal_generator_function["matrices"] = {key: np.array(new_function["matrices"][key])
                                                             for key in new_function["matrices"]}

        clear_cache()

    def insert_rates(self, new_rates):
        if len(new_rates) != 14:
            raise Exception(
                "insert_rates currently only supports the six state two binding site promoter model with 14 rate constants.")
        k_01, k_03, k_10, k_12, k_14, k_21, k_25, k_30, k_34, k_41, k_43, k_45, k_52, k_54 = new_rates
        matrices = {"c":
                        [[-(k_01), k_01, 0, 0, 0, 0],
                         [0, -(k_12), k_12, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0],
                         [0, 0, 0, -(k_34), k_34, 0],
                         [0, 0, 0, 0, -(k_45), k_45],
                         [0, 0, 0, 0, 0, 0]],
                    "1":
                        [[-(k_03), 0, 0, k_03, 0, 0],
                         [k_10, -(k_10 + k_14), 0, 0, k_14, 0],
                         [0, k_21, -(k_21 + k_25), 0, 0, k_25],
                         [k_30, 0, 0, -(k_30), 0, 0],
                         [0, k_41, 0, k_43, -(k_41 + k_43), 0],
                         [0, 0, k_52, 0, k_54, -(k_52 + k_54)]]
                    }

        new_function = {"matrices": matrices}
        self.update_infinitesimal_generator_function(new_function=new_function)

        # def call(self, in_val_dict: dict, sim_settings: dict, *args, **kwargs) -> dict:
        #     if in_val_dict != self._cur_in_val_dict:
        #         self.propensity_matrix = self._evaluate_propensity_matrix(in_val_dict)
        #         self.distribution = self._evaluate_distribution()
        #         self.entropy_production_rate = self.entropy_production_rate()
        #         self._cur_in_val_dict = in_val_dict
        #
        #     statistics = {"distribution": self.distribution,
        #                   "propensity_matrix": self.propensity_matrix,
        #                   "entropy_production_rate": self.entropy_production_rate}
        #
        #     return statistics

        # @cache_this   # Caching is not supported as it cannot handle dictionary inputs

    def __call__(self, in_val_dict, sim_settings, *args, **kwargs):
        if self._cur_in_val_dict is None or np.any([True
                                                    if len(in_val_dict[key]) != len(self._cur_in_val_dict[key])
                                                    else in_val_dict[key] != self._cur_in_val_dict[key]
                                                    for key in in_val_dict]):
            self.propensity_matrix = self._evaluate_propensity_matrix(in_val_dict, sim_settings)
            self.distribution = self._evaluate_distribution()
            self.entropy_production_rate = self._entropy_production_rate()
            self._cur_in_val_dict = in_val_dict

        statistics = {"distribution": self.distribution,
                      "propensity_matrix": self.propensity_matrix,
                      "entropy_production_rate": self.entropy_production_rate}

        return statistics

    # @cache_this
    def _evaluate_propensity_matrix(self, in_val_dict, sim_settings):
        # ToDo Adapt to linear combination of matrices scheme for efficiency reasons.
        # raise Exception("Implement Linear Combination Scheme for Propensity Matrix")

        n_samples = sim_settings["n_samples_simulation"]

        # expressions = self.infinitesimal_generator_function["expressions"]
        matrices = self.infinitesimal_generator_function["matrices"]

        propensity_mat = np.zeros((self.N_states, self.N_states, n_samples))

        for expression in matrices:
            matrix = np.expand_dims(matrices[expression], -1)
            factor_ids = expression.split("*")
            scalar = np.ones(n_samples)
            for factor_id in factor_ids:
                scalar *= in_val_dict[factor_id] if factor_id != "1" else 1

            propensity_mat += matrix * scalar

        propensity_mat = propensity_mat.transpose((2, 0, 1))
        if np.any(np.sum(propensity_mat, axis=-1) > 10 ** (-8)):
            raise Exception(
                f"Propensity Matrix is not a stochastic matrix. Row sum Zero condition is not satisfied!\nin_val_dict: {in_val_dict}\nmatrices: {matrices}\npropensity_mat: {propensity_mat}\nRow Sums: {np.sum(propensity_mat, axis=1)}")
        return propensity_mat

    # @cache_this
    def _evaluate_distribution(self):
        # distribution = np.zero(self.N_states)

        K = self.propensity_matrix
        if self.N_states == 4:
            k10, k12 = K[1, 0], K[1, 2]
            k21, k23 = K[2, 1], K[2, 3]
            k32, k30 = K[3, 2], K[3, 0]
            k01, k03 = K[0, 1], K[0, 3]

            state_weights = [k10 * k21 * k30 + k10 * k23 * k30 + k12 * k23 * k30 + k10 * k21 * k32,
                             k01 * k21 * k30 + k01 * k23 * k30 + k01 * k21 * k32 + k03 * k21 * k32,
                             k01 * k12 * k30 + k03 * k10 * k32 + k01 * k12 * k32 + k03 * k12 * k32,
                             k03 * k10 * k21 + k03 * k10 * k23 + k01 * k12 * k23 + k03 * k12 * k23]
        else:
            # ToDo Add analytical equations

            # The Steady State distribution is derived by
            #   1. Getting the propensity matrix
            #   2. Derive the vector for the left null space (apply scp.linalg.null_space() to the transposed matrix)
            #   3. Normalize the derived vector to satisfy that all probabilities add up to one

            propensity_matrix_transposed = K.transpose((0, 2, 1))

            eigvals, eigvecs = np.linalg.eig(propensity_matrix_transposed)
            # Due to numerical imprecision it happens that the eigenvalue which should be zero is not actually zero.
            eig_index = np.argmin(np.abs(eigvals), axis=-1)

            # ToDo Check the subsequent two equations
            eigvecs_transposed = eigvecs.transpose((0, 2, 1))
            state_weights = eigvecs_transposed[np.arange(eig_index.shape[0]), eig_index]

            # state_weights = state_weights.squeeze(-2)
            # We have to take the absolute values as numerical precision can potentially be insufficient and yield negative weights.
            state_weights = np.abs(state_weights)

        distribution = state_weights / np.expand_dims(np.sum(state_weights, -1), -1)
        return distribution

    # @cache_this
    def _entropy_production_rate(self):
        # Four State Code
        # if self._N_states != 4:
        #     raise Exception("Entropy Production Rate only implemented for the four state promoter")
        #
        # promoter_distribution = self.distribution(external_concentrations)
        # K = self.propensity_matrix(external_concentrations)
        # J_plus = promoter_distribution[0] * K[0, 1]
        # J_minus = promoter_distribution[1] * K[1, 0, external_concentrations]
        #
        # entropy_production_rate = (J_plus - J_minus) * np.log(J_plus / J_minus)
        # return entropy_production_rate

        distribution = self.distribution
        propensity_mat = self.propensity_matrix

        n_samples = distribution.shape[0]

        # ToDo Check correctness of this matrice!
        # flux_mat = (distribution[0] * propensity_mat[0].transpose()).transpose()
        flux_mat = propensity_mat * np.expand_dims(distribution, -1)

        upper_indices = np.triu_indices_from(flux_mat[0])
        lower_indices = np.tril_indices_from(flux_mat[0])
        # upper_vals = flux_mat[np.arange(n_samples), upper_indices]
        # lower_vals = flux_mat[np.arange(n_samples), lower_indices]
        # Prevent division by zero
        # mask = lower_vals != 0
        # upper_vals = upper_vals[mask]
        # lower_vals = lower_vals[mask]

        upper_vals = np.triu(flux_mat, k=1)
        lower_vals = np.tril(flux_mat, k=-1)
        lower_vals = lower_vals.transpose((0, 2, 1))
        mask = lower_vals != 0
        # (1 - np.tri(6, dtype=int)).flatten().repeat(n_samples).reshape(-1, n_samples).transpose().reshape(n_samples, 6,                                                                                                          -1)
        upper_vals = upper_vals[mask].reshape(n_samples, -1)
        lower_vals = lower_vals[mask].reshape(n_samples, -1)
        # np.tile(1 - np.tri(6, dtype=int), 3).transpose().reshape(3, 6, -1)

        entropy_production_rate_contributions = (upper_vals - lower_vals) * np.log(upper_vals / lower_vals)
        entropy_production_rate = np.sum(entropy_production_rate_contributions, axis=-1)

        # Reimplement this in a more efficient way

        # entropy_production_rate = 0
        #
        # for iR in range(self.N_states):
        #     for iC in range(self.N_states):
        #         if propensity_mat[0, iR, iC] == 0 or propensity_mat[0, iC, iR] == 0:
        #             continue
        #         if iC >= iR:
        #             continue
        #         flux_diff = flux_mat[0, iR, iC] - flux_mat[0, iC, iR]
        #         log_flux_diff = np.log(flux_mat[0, iR, iC] / flux_mat[0, iC, iR])
        #         entropy_production_rate += flux_diff * log_flux_diff
        #         pass

        # entropy_production_rate *= 0.5
        return entropy_production_rate

    def relaxation_time(self):
        propensity_mat = self.propensity_matrix
        eig_vals = np.linalg.eigvals(propensity_mat)

        eig_vals = np.abs(eig_vals)
        np.sort(eig_vals)
        gamma_prime = eig_vals[-1] - eig_vals[-2]  # spectral gap is equal to largest minus second largest eigenvalue
        relaxation_time = 1 / gamma_prime
        return relaxation_time

    def mixing_time_bounds(self, eps):
        promoter_distribution = self.distribution
        pi_min = np.min(promoter_distribution)
        relaxation_time = self.relaxation_time()
        lower_bound = (relaxation_time - 1) * np.log(1 / (2 * eps))
        upper_bound = relaxation_time * np.log(1 / (eps * pi_min))
        return lower_bound, upper_bound
