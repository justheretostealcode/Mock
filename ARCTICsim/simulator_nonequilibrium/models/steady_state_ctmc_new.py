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
        if in_val_dict != self._cur_in_val_dict:
            self.propensity_matrix = self._evaluate_propensity_matrix(in_val_dict)
            self.distribution = self._evaluate_distribution()
            self.entropy_production_rate = self._entropy_production_rate()
            self._cur_in_val_dict = in_val_dict

        statistics = {"distribution": self.distribution,
                      "propensity_matrix": self.propensity_matrix,
                      "entropy_production_rate": self.entropy_production_rate}

        return statistics


    # @cache_this
    def _evaluate_propensity_matrix(self, in_val_dict):
        # ToDo Adapt to linear combination of matrices scheme for efficiency reasons.
        # raise Exception("Implement Linear Combination Scheme for Propensity Matrix")

        # expressions = self.infinitesimal_generator_function["expressions"]
        matrices = self.infinitesimal_generator_function["matrices"]

        propensity_mat = np.zeros((self.N_states, self.N_states))
        for expression in matrices:
            matrix = matrices[expression]
            factor_ids = expression.split("*")
            scalar = 1
            for factor_id in factor_ids:
                scalar *= in_val_dict[factor_id] if factor_id != "1" else 1
            propensity_mat += matrix * scalar

        if any(np.sum(propensity_mat, axis=1) > 10 ** (-8)):
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

            propensity_matrix_transposed = np.transpose(K)

            eigvals, eigvecs = np.linalg.eig(propensity_matrix_transposed)
            # Due to numerical imprecision it happens that the eigenvalue which should be zero is not actually zero.
            eig_index = np.argmin(np.abs(eigvals))
            state_weights = eigvecs[:, eig_index]

            # We have to take the absolute values as numerical precision can potentially be insufficient and yield negative weights.
            state_weights = np.abs(state_weights)
        distribution = state_weights / np.sum(state_weights)

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

        # ToDo Check correctness of this matrice!
        flux_mat = (distribution * propensity_mat.transpose()).transpose()

        # Reimplement this in a more efficient way

        entropy_production_rate = 0
        # Can be simplified to n * (n-1)/2 calculations. The Factor 0.5 needs to be corrected also afterwards.
        for iR in range(self.N_states):
            for iC in range(self.N_states):
                if propensity_mat[iR, iC] == 0 or propensity_mat[iC, iR] == 0:
                    continue
                if iC >= iR:
                    continue
                flux_diff = flux_mat[iR, iC] - flux_mat[iC, iR]
                log_flux_diff = np.log(flux_mat[iR, iC] / flux_mat[iC, iR])
                entropy_production_rate += flux_diff * log_flux_diff
                pass

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
