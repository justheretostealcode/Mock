"""
Author: Erik Kubaczka
E-Mail: erik.kubaczka@tu-darmstadt.de
"""

import numpy as np

from ARCTICsim.simulator_nonequilibrium.models.steady_state_ctmc_new import SteadyStateCTMC


class PromoterModel(SteadyStateCTMC):

    def __init__(self,
                 num_states,
                 num_trainable_parameters,
                 infinitesimal_generator_function,
                 per_state_promoter_activity):
        super().__init__(num_states, infinitesimal_generator_function)

        self.num_trainable_parameters = num_trainable_parameters

        self._promoter_activity = np.array(per_state_promoter_activity)

    def __call__(self, in_val_dict: dict, sim_settings: dict, *args, **kwargs) -> dict:
        """
        :param in_val:
        :param args:
        :param kwargs:
        :return:
        """

        statistics = super()(in_val_dict, sim_settings=sim_settings)

        results = dict(statistics)

        results["average_promoter_activity"] = self._promoter_activity @ statistics["distribution"]
        results["average_promoter_activity_per_state"] = list(self._promoter_activity * statistics["distribution"])
        results["energy_dissipation_rate"] = statistics["entropy_production_rate"]
        results["energy_p"] = statistics["entropy_production_rate"]
        results["promoter_activity_per_state"] = self._promoter_activity

        return results
