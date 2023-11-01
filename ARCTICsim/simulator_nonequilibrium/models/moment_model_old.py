"""
Author: Erik Kubaczka
E-Mail: erik.kubaczka@tu-darmstadt.de
"""

import numpy as np


class Moments:

    def __init__(self,
                 N_promoter_states,
                 k_promoter,
                 transcription_rates,
                 mRNA_degradation_rate,
                 translation_rate,
                 protein_degradation_rate,
                 normalize_cycle_time=False,
                 c=-1):
        def str_to_expr(c, r0, l12, l21, l23, l32, l34, l43, l41, l14, str_expr):
            r_0 = r0
            val = eval(str_expr.replace("^", "**"))
            return val

        if N_promoter_states != 4 or len(transcription_rates) != 4:
            raise Exception("The current implementation of moments supports only four state promoter models.")

        self._N_promoter_states = N_promoter_states
        self._K = k_promoter  # k_promoter depends on the ligand concentration c
        self._mu1_rates = transcription_rates
        self._d1 = mRNA_degradation_rate
        self._mu2 = translation_rate
        self._d2 = protein_degradation_rate

        if normalize_cycle_time and c != -1:
            tau_b = self.cycle_time(c=c)
            self._K = lambda i1, i2, c: tau_b * k_promoter(i1, i2, c)

        with open("symbolic_expressions/sharpness_expression.txt", "r") as file:
            str_lammers_sharpness = file.readline()
        self.lammers_sharpness_equation = lambda c, r0, l12, l21, l23, l32, l34, l43, l41, l14: str_to_expr(c, r0, l12,
                                                                                                            l21, l23,
                                                                                                            l32, l34,
                                                                                                            l43, l41,
                                                                                                            l14,
                                                                                                            str_lammers_sharpness)

    @property
    def K(self):
        return self._K

    @K.setter
    def K(self, value):
        self._K = value

    @property
    def transcription_rates(self):
        return self._mu1_rates

    @transcription_rates.setter
    def transcription_rates(self, value):
        self._mu1_rates = value

    @property
    def mRNA_degradation_rate(self):
        return self._d1

    @mRNA_degradation_rate.setter
    def mRNA_degradation_rate(self, value):
        self._d1 = value

    @property
    def translation_rate(self):
        return self._mu2

    @translation_rate.setter
    def translation_rate(self, value):
        self._mu2 = value

    @property
    def protein_degradation_rate(self):
        return self._d2

    @protein_degradation_rate.setter
    def protein_degradation_rate(self, value):
        self._d2 = value

    def promoter_distribution(self, c):
        state_weights = [self._K(1, 0, c) * self._K(2, 1, c) * self._K(3, 0, c)
                         + self._K(1, 0, c) * self._K(2, 3, c) * self._K(3, 0, c)
                         + self._K(1, 2, c) * self._K(2, 3, c) * self._K(3, 0, c)
                         + self._K(1, 0, c) * self._K(2, 1, c) * self._K(3, 2, c),
                         self._K(0, 1, c) * self._K(2, 1, c) * self._K(3, 0, c)
                         + self._K(0, 1, c) * self._K(2, 3, c) * self._K(3, 0, c)
                         + self._K(0, 1, c) * self._K(2, 1, c) * self._K(3, 2, c)
                         + self._K(0, 3, c) * self._K(2, 1, c) * self._K(3, 2, c),
                         self._K(0, 1, c) * self._K(1, 2, c) * self._K(3, 0, c)
                         + self._K(0, 3, c) * self._K(1, 0, c) * self._K(3, 2, c)
                         + self._K(0, 1, c) * self._K(1, 2, c) * self._K(3, 2, c)
                         + self._K(0, 3, c) * self._K(1, 2, c) * self._K(3, 2, c),
                         self._K(0, 3, c) * self._K(1, 0, c) * self._K(2, 1, c)
                         + self._K(0, 3, c) * self._K(1, 0, c) * self._K(2, 3, c)
                         + self._K(0, 1, c) * self._K(1, 2, c) * self._K(2, 3, c)
                         + self._K(0, 3, c) * self._K(1, 2, c) * self._K(2, 3, c)]

        promoter_distribution = state_weights / np.sum(state_weights)

        return promoter_distribution

    def propensity_mat(self, c):
        propensity_mat = np.array(
            [[self._K(iR, iC, c) for iC in range(self._N_promoter_states)] for iR in range(self._N_promoter_states)])
        return propensity_mat

    def cycle_time(self, c):
        a = self._mu1_rates > 0
        off_state_counts = np.count_nonzero(1 - a)
        on_state_counts = np.count_nonzero(a)

        promoter_distribution = self.promoter_distribution(c)

        propensity_mat = self.propensity_mat(c)

        K_off_mat = np.zeros(shape=(off_state_counts + 1, off_state_counts + 1))
        for iR in range(off_state_counts):
            for iC in range(off_state_counts):
                if iR == iC:
                    continue
                K_off_mat[iR, iC] = propensity_mat[iR, iC]

        for iI in range(off_state_counts):
            K_off_mat[off_state_counts, iI] = np.sum([a[iR] * propensity_mat[iR, iI]
                                                      for iR in range(self._N_promoter_states)])
            K_off_mat[iI, off_state_counts] = np.sum([a[iR] * propensity_mat[iI, iR]
                                                      for iR in range(self._N_promoter_states)])

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
                                          for iR in range(self._N_promoter_states)])
            K_on_mat[1 + iI, 0] = np.sum([(1 - a[iR]) * propensity_mat[off_state_counts + iI, iR]
                                          for iR in range(self._N_promoter_states)])

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
        sharpness = self.lammers_sharpness_equation(c=c, r0=1, l12=propensity_mat[0, 1], l21=propensity_mat[1, 0],
                                                    l23=propensity_mat[1, 2], l32=propensity_mat[2, 1],
                                                    l34=propensity_mat[2, 3], l43=propensity_mat[3, 2],
                                                    l41=propensity_mat[3, 0], l14=propensity_mat[0, 3])
        return sharpness

    def relaxation_time(self, c):
        propensity_mat = self.propensity_mat(c)
        eig_vals, eig_vecs = np.linalg.eig(propensity_mat)

        eig_vals = np.abs(eig_vals)
        np.sort(eig_vals)
        gamma_prime = eig_vals[-1] - eig_vals[-2]  # spectral gap is equal to largest minus second largest eigen value
        relaxation_time = 1 / gamma_prime
        return relaxation_time

    def mixing_time_bounds(self, c, eps):
        promoter_distribution = self.promoter_distribution(c)
        pi_min = np.min(promoter_distribution)
        relaxation_time = self.relaxation_time(c)
        lower_bound = (relaxation_time - 1) * np.log(1 / (2 * eps))
        upper_bound = relaxation_time * np.log(1 / (eps * pi_min))
        return lower_bound, upper_bound

    def entropy_production_rate(self, c):
        if self._N_promoter_states != 4:
            raise Exception("Entropy Production Rate only implemented for the four state promoter")

        promoter_distribution = self.promoter_distribution(c)
        K = self.propensity_mat(c)
        J_plus = promoter_distribution[0] * K[0, 1]
        J_minus = promoter_distribution[1] * K[1, 0, c]

        entropy_production_rate = (J_plus - J_minus) * np.log(J_plus / J_minus)
        return entropy_production_rate

    def E_mu1(self, c):
        vals = self.promoter_distribution(c) * self._mu1_rates
        E_mu1 = np.sum(vals)
        return E_mu1

    def E_M(self, c):
        E_M = self.E_mu1(c) / self._d1
        return E_M

    def E_P(self, c):
        E_M = self.E_M(c)
        E_P = E_M * self.translation_rate / self.protein_degradation_rate
        return E_P

    def E_M_Xx(self, c):
        LAMBDA = np.array([[float(self._K(iR, iC, c)) for iC in range(self._N_promoter_states)]
                           for iR in range(self._N_promoter_states)])

        M = self._d1 * np.eye(self._N_promoter_states) - np.transpose(LAMBDA)
        M_inv = np.linalg.inv(M)

        b = self._mu1_rates * self.promoter_distribution(c) * 1
        E_M_Xx = M_inv @ b

        # E_M_Xx = np.array([0.01,  0.049, 2.238, 1.139])
        return E_M_Xx

    def E_M_mu1(self, c):
        vals = self._mu1_rates * self.E_M_Xx(c)
        E_M_mu1 = np.sum(vals)
        return E_M_mu1

    def E_P_Xx(self, c):
        LAMBDA = np.array([[float(self._K(iR, iC, c)) for iC in range(self._N_promoter_states)]
                           for iR in range(self._N_promoter_states)])

        Mat = self._d2 * np.eye(self._N_promoter_states) - np.transpose(LAMBDA)
        M_inv = np.linalg.inv(Mat)

        E_M_Xx = self.E_M_Xx(c)

        b = E_M_Xx * 1
        E_P_Xx = self._mu2 * M_inv @ b

        return E_P_Xx

    def E_P_mu1(self, c):
        E_P_Xx = self.E_P_Xx(c)
        vals = self._mu1_rates * E_P_Xx
        E_P_mu1 = np.sum(vals)

        return E_P_mu1

    def E_M_P(self, c):
        E_P_mu1 = self.E_P_mu1(c)
        E_M_squared = self.E_M_squared(c)

        E_M_P = E_P_mu1 + self._mu2 * E_M_squared
        E_M_P = E_M_P / (self._d1 + self._d2)
        return E_M_P

    def E_M_squared(self, c):
        E_mu1 = self.E_mu1(c)
        E_M_mu1 = self.E_M_mu1(c)
        E_M = self.E_M(c)

        E_M_squared = E_mu1 + 2 * E_M_mu1 + self._d1 * E_M
        E_M_squared = E_M_squared / (2 * self._d1)

        return E_M_squared

    def E_P_squared(self, c):
        E_M = self.E_M(c)
        E_M_P = self.E_M_P(c)
        E_P = self.E_P(c)

        E_P_Squared = self._mu2 * E_M + 2 * self._mu2 * E_M_P + self._d2 * E_P
        E_P_Squared = E_P_Squared / (2 * self._d2)

        return E_P_Squared

    def mean_M(self, c):
        return self.E_M(c)

    def var_M(self, c):
        return self.E_M_squared(c) - (self.E_M(c)) ** 2

    def mean_P(self, c):
        return self.E_P(c)

    def var_P(self, c):
        return self.E_P_squared(c) - (self.E_P(c)) ** 2



"""
Helper Functions
"""



def rates_to_lammers_rates(rates):
    k12, k21, k23, k32, k34, k43, k41, k14 = rates

    k_i = k41
    k_a = k14
    k_b = k12
    k_u = k21

    n_ab = k23 / k14
    n_ib = k32 / k41
    n_ua = k34 / k21
    n_ba = k43 / k12

    return k_a, k_i, k_b, k_u, n_ba, n_ib, n_ua, n_ab


def lammers_rates_to_rates(lammers_rates):
    k_a, k_i, k_b, k_u, n_ba, n_ib, n_ua, n_ab = lammers_rates

    k12 = k_b
    k21 = k_u

    k23 = n_ab * k_a
    k32 = n_ib * k_i

    k34 = n_ua * k_u
    k43 = n_ba * k_b

    k41 = k_i
    k14 = k_a

    return k12, k21, k23, k32, k34, k43, k41, k14


def get_propensity_function(k_b, k_u, k_a, k_i, n_ab, n_ua, n_ba, n_ib, atp=1, adp_pi=1, c_func=None):
    if c_func is None:
        c_func = lambda l: l
        # gamma = 100
        # c_func = lambda l: 1 + gamma * l

    k12 = lambda l: c_func(l) * k_b
    k21 = lambda l: k_u
    k23 = lambda l: n_ab * k_a
    k32 = lambda l: n_ib * k_i
    k34 = lambda l: n_ua * k_u
    k43 = lambda l: c_func(l) * n_ba * k_b
    k41 = lambda l: atp * k_i
    k14 = lambda l: adp_pi * k_a

    k11 = lambda l: - c_func(l) * k_b - adp_pi * k_a
    k22 = lambda l: - k_u - n_ab * k_a
    k33 = lambda l: - n_ib * k_i - n_ua * k_u
    k44 = lambda l: - atp * k_i - c_func(l) * n_ba * k_b

    zero_func = lambda l: 0

    promoter_rates = [[k11, k12, zero_func, k14],
                      [k21, k22, k23, zero_func],
                      [zero_func, k32, k33, k34],
                      [k41, zero_func, k43, k44]]

    K = lambda i1, i2, l: promoter_rates[i1][i2](l)

    return K
