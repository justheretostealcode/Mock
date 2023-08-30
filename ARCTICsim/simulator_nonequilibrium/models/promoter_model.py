"""
Author: Erik Kubaczka
E-Mail: erik.kubaczka@tu-darmstadt.de
"""


import numpy as np

from models.steady_state_ctmc import SteadyStateCTMC

from ARCTICsim.simulator_nonequilibrium.models.custom_cache import cache_this


class PromoterModel(SteadyStateCTMC):
    _LAMMERS_SHARPNESS_EXPRESSION = "c*((r_0*(l12*l23*l41 + l14*l21*l43 + l14*l23*l43 + 2*c*l12*l23*l43))/(l14*l21*l32 + l14*l21*l34 + l14*l23*l34 + l21*l32*l41 + l21*l34*l41 + l23*l34*l41 + c*l12*l23*l34 + c*l12*l23*l41 + c*l14*l21*l43 + c*l14*l23*l43 + c*l12*l32*l41 + c*l12*l34*l41 + c*l14*l32*l43 + c*l21*l32*l43 + c^2*l12*l23*l43 + c^2*l12*l32*l43) - (r_0*(l14*l21*l32 + l14*l21*l34 + l14*l23*l34 + c*l12*l23*l34)*(l12*l23*l34 + l12*l23*l41 + l14*l21*l43 + l14*l23*l43 + l12*l32*l41 + l12*l34*l41 + l14*l32*l43 + l21*l32*l43 + 2*c*l12*l23*l43 + 2*c*l12*l32*l43))/(l14*l21*l32 + l14*l21*l34 + l14*l23*l34 + l21*l32*l41 + l21*l34*l41 + l23*l34*l41 + c*l12*l23*l34 + c*l12*l23*l41 + c*l14*l21*l43 + c*l14*l23*l43 + c*l12*l32*l41 + c*l12*l34*l41 + c*l14*l32*l43 + c*l21*l32*l43 + c^2*l12*l23*l43 + c^2*l12*l32*l43)^2 - (r_0*(c*l12*l23*l41 + c*l14*l21*l43 + c*l14*l23*l43 + c^2*l12*l23*l43)*(l12*l23*l34 + l12*l23*l41 + l14*l21*l43 + l14*l23*l43 + l12*l32*l41 + l12*l34*l41 + l14*l32*l43 + l21*l32*l43 + 2*c*l12*l23*l43 + 2*c*l12*l32*l43))/(l14*l21*l32 + l14*l21*l34 + l14*l23*l34 + l21*l32*l41 + l21*l34*l41 + l23*l34*l41 + c*l12*l23*l34 + c*l12*l23*l41 + c*l14*l21*l43 + c*l14*l23*l43 + c*l12*l32*l41 + c*l12*l34*l41 + c*l14*l32*l43 + c*l21*l32*l43 + c^2*l12*l23*l43 + c^2*l12*l32*l43)^2 + (l12*l23*l34*r_0)/(l14*l21*l32 + l14*l21*l34 + l14*l23*l34 + l21*l32*l41 + l21*l34*l41 + l23*l34*l41 + c*l12*l23*l34 + c*l12*l23*l41 + c*l14*l21*l43 + c*l14*l23*l43 + c*l12*l32*l41 + c*l12*l34*l41 + c*l14*l32*l43 + c*l21*l32*l43 + c^2*l12*l23*l43 + c^2*l12*l32*l43))"

    def __init__(self,
                 num_states,
                 num_trainable_parameters,
                 infinitesimal_generator_function,
                 transcription_rates,
                 normalize_cycle_time=False,
                 external_concentrations_for_normalization=-1):
        super().__init__(num_states, infinitesimal_generator_function)

        self.num_trainable_parameters = num_trainable_parameters

        def str_to_expr(c, r0, l12, l21, l23, l32, l34, l43, l41, l14, str_expr):
            r_0 = r0
            val = eval(str_expr.replace("^", "**"))
            return val

        self._mu1_rates = transcription_rates

        if normalize_cycle_time and external_concentrations_for_normalization >= 0:
            tau_b = self.cycle_time(external_concentrations=external_concentrations_for_normalization)
            super().infinitesimal_generator_function = lambda i1, i2, external_concentrations: \
                tau_b * super().infinitesimal_generator_function(i1, i2, external_concentrations)

        self.lammers_sharpness_equation = lambda c, r0, l12, l21, l23, l32, l34, l43, l41, l14: \
            str_to_expr(c, r0, l12, l21, l23, l32, l34, l43, l41, l14, self._LAMMERS_SHARPNESS_EXPRESSION)
        pass

    @cache_this
    def E_mu1(self, external_concentrations):
        vals = super().distribution(external_concentrations) * self._mu1_rates
        E_mu1 = np.sum(vals)
        return E_mu1

    def E_mu1x_Xx_squared(self, external_concentrations):
        E_mu1_X_squared = self._mu1_rates * self.distribution(external_concentrations) * 1
        return E_mu1_X_squared

    def cycle_time(self, external_concentrations):
        # ToDo Refactor this method to remove duplicate code (e.g. follow DRY rule)
        # Approach for ET_OFF_ON and ET_ON_OFF is equal except the states considered
        # One could utilize the SteadyStateCTMC class for the steady state distributions etc.
        a = self._mu1_rates > 0
        off_state_counts = np.count_nonzero(1 - a)
        on_state_counts = np.count_nonzero(a)

        promoter_distribution = super().distribution(external_concentrations)

        propensity_mat = self.propensity_matrix(external_concentrations)

        K_off_mat = np.zeros(shape=(off_state_counts + 1, off_state_counts + 1))
        for iR in range(off_state_counts):
            for iC in range(off_state_counts):
                if iR == iC:
                    continue
                K_off_mat[iR, iC] = propensity_mat[iR, iC]

        for iI in range(off_state_counts):
            K_off_mat[off_state_counts, iI] = np.sum([a[iR] * propensity_mat[iR, iI]
                                                      for iR in range(super().N_states)])
            K_off_mat[iI, off_state_counts] = np.sum([a[iR] * propensity_mat[iI, iR]
                                                      for iR in range(super().N_states)])

        if any(np.diag(K_off_mat) != 0):
            print("Diagonal should be 0")

        diag_vals = np.sum(K_off_mat, axis=1)
        K_off_mat = K_off_mat - np.diag(diag_vals)

        # for iD in range(off_state_counts + 1):
        #    K_off_mat[iD, iD] = - np.sum([K_off_mat[iR, iD] for iR in range(self._N_promoter_states)])
        if not all(np.abs(np.sum(K_off_mat, axis=1)) <= 10 ** (-10)):
            # print("Unsuccesful creation of K_off_mat")
            pass

        K_on_mat = np.zeros(shape=(on_state_counts + 1, on_state_counts + 1))
        for iR in range(off_state_counts):
            for iC in range(off_state_counts):
                if iR == iC:
                    continue
                K_on_mat[1 + iR, 1 + iC] = propensity_mat[off_state_counts + iR, off_state_counts + iC]

        for iI in range(on_state_counts):
            K_on_mat[0, 1 + iI] = np.sum([(1 - a[iR]) * propensity_mat[iR, off_state_counts + iI]
                                          for iR in range(super().N_states)])
            K_on_mat[1 + iI, 0] = np.sum([(1 - a[iR]) * propensity_mat[off_state_counts + iI, iR]
                                          for iR in range(super().N_states)])

        if any(np.diag(K_on_mat) != 0):
            print("Diagonal should be 0")

        diag_vals = np.sum(K_on_mat, axis=1)
        K_on_mat = K_on_mat - np.diag(diag_vals)

        # for iD in range(on_state_counts + 1):
        #     K_on_mat[iD, iD] = - np.sum([K_on_mat[iR, iD] for iR in range(self._N_promoter_states)])

        if not all(np.abs(np.sum(K_on_mat, axis=1)) <= 10 ** (-8)):
            # print("Unsuccesful creation of K_on_mat")
            pass

        ####
        # Procedure for ET_OFF_ON
        ####
        eigvals, eigvecs = np.linalg.eig(np.transpose(K_off_mat))
        eigvals = np.abs(eigvals)
        eigvals = eigvals / np.sum(eigvals)  # Normalize to ensure that the thresholding for 0 values works properly
        zero_eig_index = np.where(eigvals < 10 ** (-12))[0]
        if len(zero_eig_index) == 0:
            print("There is no zero eig val")

        distribution_off = eigvecs[:, zero_eig_index[0]]
        distribution_off = np.abs(distribution_off)  # Used to prevent imaginery part
        distribution_off = distribution_off / np.sum(distribution_off)

        PI = np.repeat(distribution_off, off_state_counts + 1).reshape(off_state_counts + 1, -1).transpose()

        Z = np.linalg.inv(PI - K_off_mat) - PI

        # Equation S20 of Lammers
        et_OFF_to_ON = (Z[off_state_counts, off_state_counts] - Z[:off_state_counts, off_state_counts]) / \
                       distribution_off[off_state_counts]

        f_OFF = propensity_mat[:off_state_counts] @ (a * promoter_distribution)

        ET_OFF_ON = np.sum(et_OFF_to_ON * f_OFF) / np.sum(f_OFF)

        ####
        # Procedure for ET_ON_OFF
        ####
        eigvals, eigvecs = np.linalg.eig(np.transpose(K_on_mat))
        eigvals = np.abs(eigvals)
        eigvals = eigvals / np.sum(eigvals)  # Normalize to ensure that the thresholding for 0 values works properly
        zero_eig_index = np.where(eigvals < 10 ** (-12))[0]
        if len(zero_eig_index) == 0:
            print("There is no zero eig val")

        distribution_on = eigvecs[:, zero_eig_index[0]]
        distribution_on = np.abs(distribution_on)  # Used to prevent imaginery part
        distribution_on = distribution_on / np.sum(distribution_on)

        PI = np.repeat(distribution_on, on_state_counts + 1).reshape(on_state_counts + 1, -1).transpose()

        # This is the correct formular for row sum zero matrix
        Z = np.linalg.inv(PI - K_on_mat) - PI

        # Equation S20 of Lammers2023 (Unclear whether this matches their definition as they use a column sum zero matrix (column to row transitions))
        et_ON_to_OFF = (Z[0, 0] - Z[1:, 0]) / distribution_on[1:]

        f_ON = propensity_mat[:off_state_counts] @ (a * promoter_distribution)

        ET_ON_OFF = np.sum(et_ON_to_OFF * f_ON) / np.sum(f_ON)

        burst_cycle_time = ET_OFF_ON + ET_ON_OFF

        return burst_cycle_time

    # Correct version of Lammers precision accounting for row sum zero propensity matrix
    def lammers_precision(self, c):
        a = self._mu1_rates > 0

        promoter_distribution = self.promoter_distribution(c)
        propensity_mat = self.propensity_mat(c)

        PI = np.repeat(promoter_distribution, self._N_promoter_states).reshape(self._N_promoter_states, -1).transpose()

        Z = np.linalg.inv(PI - propensity_mat) - PI

        var = 2 * (a * promoter_distribution) @ Z @ a

        std_dev = np.sqrt(var)
        return 1 / std_dev

    def lammers_sharpness(self, c):
        propensity_mat = self.propensity_mat(c)
        # This equation already includes the multiplication by c
        sharpness = self.lammers_sharpness_equation(c=c, r0=1,
                                                    l12=propensity_mat[0, 1], l21=propensity_mat[1, 0],
                                                    l23=propensity_mat[1, 2], l32=propensity_mat[2, 1],
                                                    l34=propensity_mat[2, 3], l43=propensity_mat[3, 2],
                                                    l41=propensity_mat[3, 0], l14=propensity_mat[0, 3])
        return sharpness
