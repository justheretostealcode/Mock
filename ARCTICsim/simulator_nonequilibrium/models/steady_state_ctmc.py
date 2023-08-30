"""
Author: Erik Kubaczka
E-Mail: erik.kubaczka@tu-darmstadt.de
"""

import itertools

import numpy as np

from ARCTICsim.simulator_nonequilibrium.models.custom_cache import cache_this, clear_cache


class SteadyStateCTMC:

    def __init__(self,
                 N_states,
                 infinitesimal_generator_function):
        self.N_states = N_states
        self.update_infinitesimal_generator_function(infinitesimal_generator_function)
        # self.infinitesimal_generator_function = infinitesimal_generator_function

        pass

    def update_infinitesimal_generator_function(self, new_function):
        self.infinitesimal_generator_function = new_function

        clear_cache()
        pass

    @cache_this
    def propensity_matrix(self, external_concentrations):

        # propensity_mat_legacy = np.array([[self.infinitesimal_generator_function(iR, iC, external_concentrations)
        #                             for iC in range(self.N_states)]
        #                            for iR in range(self.N_states)], dtype=float)
        propensity_mat = np.array([self.infinitesimal_generator_function(index[0], index[1], external_concentrations)
                                   for index in itertools.product(range(self.N_states), range(self.N_states))],
                                  dtype=float)
        propensity_mat = propensity_mat.reshape(self.N_states, -1)
        return propensity_mat

    @cache_this
    def distribution(self, external_concentrations):

        # The Steady State distribution is derived by
        #   1. Getting the propensity matrix
        #   2. Derive the vector for the left null space (apply scp.linalg.null_space() to the transposed matrix)
        #   3. Normalize the derived vector to satisfy that all probabilities add up to one

        propensity_matrix = self.propensity_matrix(external_concentrations=external_concentrations)

        propensity_matrix_transposed = np.transpose(propensity_matrix)

        eigvals, eigvecs = np.linalg.eig(propensity_matrix_transposed)
        # Thresholding is required due to successfully identify
        eps_float = 2.220446049250313e-16
        # zero_threshold = eps_float * np.product(propensity_matrix_transposed.shape) * np.max(
        #     propensity_matrix_transposed)
        # eig_index = np.where(np.abs(eigvals) < zero_threshold)[0][0]
        eig_index = np.argmin(np.abs(eigvals))
        state_weights = eigvecs[:, eig_index]

        # the Null Space method does not work as it happens to produce negative probabilities
        # the function null_space
        # considers all singular values $s$ to be zero in case they are smaller than rcond * max(s)
        # null_space_basis = scp.linalg.null_space(propensity_matrix_transposed, rcond=None)
        # if null_space_basis.shape[1] == 0:
        #     raise Exception("Unable to derive steady state distribution as null space of matrix is empty.")
        #
        # state_weights = null_space_basis[:, 0]

        # Currently required as the eigenvalue method is of approximative nature only
        state_weights = np.abs(state_weights)
        distribution = state_weights / np.sum(state_weights)

        return distribution

    def relaxation_time(self, external_concentrations):
        propensity_mat = self.propensity_matrix(external_concentrations=external_concentrations)
        eig_vals = np.linalg.eigvals(propensity_mat)
        # The eigenvalues of

        eig_vals = np.abs(eig_vals)
        np.sort(eig_vals)
        gamma_prime = eig_vals[-1] - eig_vals[-2]  # spectral gap is equal to largest minus second largest eigen value
        relaxation_time = 1 / gamma_prime
        return relaxation_time

    def mixing_time_bounds(self, external_concentrations, eps):
        promoter_distribution = self.distribution(external_concentrations)
        pi_min = np.min(promoter_distribution)
        relaxation_time = self.relaxation_time(external_concentrations)
        lower_bound = (relaxation_time - 1) * np.log(1 / (2 * eps))
        upper_bound = relaxation_time * np.log(1 / (eps * pi_min))
        return lower_bound, upper_bound

    @cache_this
    def entropy_production_rate(self, external_concentrations):
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

        distribution = self.distribution(external_concentrations=external_concentrations)
        propensity_mat = self.propensity_matrix(external_concentrations=external_concentrations)

        flux_mat = (distribution * propensity_mat.transpose()).transpose()

        entropy_production_rate = 0
        # Can be simplified to n * (n-1)/2 calculations. The Factor 0.5 needs to be corrected also afterwards.
        for iR in range(self.N_states):
            for iC in range(self.N_states):
                if propensity_mat[iR, iC] == 0 or propensity_mat[iC, iR] == 0:
                    continue
                flux_diff = flux_mat[iR, iC] - flux_mat[iC, iR]
                log_flux_diff = np.log10(flux_mat[iR, iC] / flux_mat[iC, iR])
                entropy_production_rate += flux_diff * log_flux_diff
                pass

        entropy_production_rate *= 0.5
        return entropy_production_rate
