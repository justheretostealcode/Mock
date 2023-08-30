"""
Author: Erik Kubaczka
E-Mail: erik.kubaczka@tu-darmstadt.de
"""
import scipy as scp
import numpy as np

from models.energy_aware_moment_model import EnergyAwareMomentModel


class FourStatePromoterModel(EnergyAwareMomentModel):
    num_trainable_parameters = 8

    def __init__(self,
                 gate_params, rates=None):

        model_params = dict(gate_params)
        promoter_model = model_params["GENE"]

        c_func_str = (promoter_model["C_FUNC"])
        c_func = eval(c_func_str)
        if "INFINITESIMAL_GENERATORS" in promoter_model:
            infinitesimal_generator_function = None
            raise Exception("Parse Infinitesimal generators")
        elif rates is not None:
            best_params = list(map(float, rates))
            model_params["GENE"]["RATES"] = FourStatePromoterModel.rates_to_rates_dict(best_params)
            infinitesimal_generator_function = self.get_infinitesimal_generator_function(rates=rates, c_func=c_func)
        elif "RATES" in promoter_model:
            rates_dict = promoter_model["RATES"]
            rates = FourStatePromoterModel.rates_dict_to_rates(rates_dict)
            infinitesimal_generator_function = self.get_infinitesimal_generator_function(rates=rates, c_func=c_func)
        else:
            # ToDo Maybe change that this is allowed in order to allow for the addition of tunable parameters afterwards
            raise Exception(
                "?? ToDo ?? You need to specify either the rates or the infinitesimal generators in the gate library")

        super().__init__(model_params,
                         num_trainable_parameters=FourStatePromoterModel.num_trainable_parameters,
                         infinitesimal_generator_function=infinitesimal_generator_function)
        # EnergyAwareMomentModel(gate_parameters,
        #                       infinitesimal_generator_function=infinitesimal_generator_function)

    """                                                                                                                                      
       Helper Functions                                                                                                                         
       """

    def rates_to_lammers_rates(self, rates):
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

    def lammers_rates_to_rates(self, lammers_rates):
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

    def get_infinitesimal_generator_function(cls, rates=None, lammers_rates=None, atp=1, adp_pi=1, c_func=None):
        def infinitesimal_generator_func(i1, i2, l):
            return promoter_rates[i1][i2](l)
        if c_func is None:
            c_func = lambda l: l
            # gamma = 100
            # c_func = lambda l: 1 + gamma * l

        if rates is None and lammers_rates is None:
            raise Exception("One of the rate variables needs to be defined")

        if lammers_rates is not None:
            cur_rates = cls.lammers_rates_to_rates(lammers_rates=rates)
        if rates is not None:
            cur_rates = rates

        k12_rc, k21_rc, k23_rc, k32_rc, k34_rc, k43_rc, k41_rc, k14_rc = cur_rates

        k12 = lambda l: c_func(l) * k12_rc
        k21 = lambda l: k21_rc
        k23 = lambda l: k23_rc
        k32 = lambda l: k32_rc
        k34 = lambda l: k34_rc
        k43 = lambda l: c_func(l) * k43_rc
        k41 = lambda l: atp * k41_rc
        k14 = lambda l: adp_pi * k14_rc

        # k11 = lambda l: - c_func(l) * k_b - adp_pi * k_a
        # k22 = lambda l: - k_u - n_ab * k_a
        # k33 = lambda l: - n_ib * k_i - n_ua * k_u
        # k44 = lambda l: - atp * k_i - c_func(l) * n_ba * k_b

        k11 = lambda l: - k12(l) - k14(l)
        k22 = lambda l: - k21(l) - k23(l)
        k33 = lambda l: - k32(l) - k34(l)
        k44 = lambda l: - k41(l) - k43(l)

        zero_func = lambda l: 0

        promoter_rates = [[k11, k12, zero_func, k14],
                          [k21, k22, k23, zero_func],
                          [zero_func, k32, k33, k34],
                          [k41, zero_func, k43, k44]]

        #K = lambda i1, i2, l: promoter_rates[i1][i2](l)
        K = infinitesimal_generator_func

        return K

    @staticmethod
    def rates_dict_to_rates(rates_dict):
        rate_ids = ["k12", "k21", "k23", "k32", "k34", "k43", "k41", "k14"]
        rates = [rates_dict[rate_id] for rate_id in rate_ids]
        rates = np.array(rates)
        if any(np.isnan(rates)):
            raise ValueError("None of the rates is allowed to be nan.")
        if any(np.isinf(rates)):
            raise ValueError("None of the rates is allowed to be infinity.")
        if any(rates < 0):
            raise ValueError("None of the rates is allowed to be negative.")
        return rates

    @staticmethod
    def rates_to_rates_dict(rates):
        rate_ids = ["k12", "k21", "k23", "k32", "k34", "k43", "k41", "k14"]
        rates_dict = {rate_id: val for rate_id, val in zip(rate_ids, rates)}
        return rates_dict

    @staticmethod
    def insert_rates_in_gate_entry(gate_entry, rates):
        gate_entry_intermediate = dict(gate_entry)
        rates_dict = FourStatePromoterModel.rates_to_rates_dict(rates)
        gate_entry_intermediate["GENE"]["RATES"] = rates_dict
        return gate_entry_intermediate


######## Class End


{
    "k12": 0.030442802032205233,
    "k21": 86.81272257909148,
    "k23": 0.16517146313599387,
    "k32": 1000.0,
    "k34": 36.185097777762266,
    "k43": 0.007813606538304319,
    "k41": 32.41336059069566,
    "k14": 3.3243041092547623

}
(1, 10, 10000, 1e-05, 0.001, 10000, 1, 10)
EXAMPLE_GATE_DATA = {
    "GENE": {
        "N_PROMOTER_STATES": 4,
        "C_FUNC": "lambda l: 1 + 10 * l",
        "RATES": {
            "k12": 1,
            "k21": 10,
            "k23": 10000,
            "k32": 1e-05,
            "k34": 0.001,
            "k43": 10000,
            "k41": 1,
            "k14": 10
        },
    },
    "mRNA": {
        "TRANSCRIPTION_RATES": [
            0.0,
            0.0,
            40.0,
            40.0
        ],
        "mRNA_DEGRADATION_RATE": 5.83e-4,
        "e": 16,
        "e_const": 0,
        "l": 752.0
    },
    "PROTEIN": {
        "TRANSLATION_RATE": 6,
        "PROTEIN_DEGRADATION_RATE": 2.83e-4,
        "e": 42,
        "e_const": 52,
        "l": 210.0
    }
}

if __name__ == '__main__':
    print("Started")
    gate = FourStatePromoterModel(EXAMPLE_GATE_DATA)
    # c = 10 ** (-10)
    c = 100
    mean_p = gate.mean_P(c)
    var_p = gate.var_P(c)
    energy_rate = gate.energy_rate(c)

    print("Mean P:", mean_p)
    print("Var P: ", var_p)
    print("Energy Rate:", energy_rate)
    pass
