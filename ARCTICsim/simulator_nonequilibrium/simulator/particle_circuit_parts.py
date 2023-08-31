"""
Author: Erik Kubaczka
"""
import json
from collections import OrderedDict

from models.four_state_promoter_model import FourStatePromoterModel
import numpy as np


class Device:
    def __init__(self, id):
        self.identifier = id
        self.energy_rate = 0

    def __str__(self):
        return f"{self.__class__.__name__} ({self.identifier})"


class InputOutput(Device):
    def __call__(self, val:float, sim_settings) -> float:
        return val


class Input(InputOutput):
    def __init__(self, id):
        super().__init__(id)
        self.type = "INPUT"


class LutInput(InputOutput):
    def __init__(self, input_entry):
        id = input_entry["identifier"]
        super().__init__(id)

        self.values = input_entry["biorep"]
        self.type = "INPUT"

        self.entry = input_entry

    def __call__(self, bool_val: int, sim_settings) -> float:
        # ToDo Remove Conversion to RPU
        output = self.values[str(bool_val)]
        out_val = output  # / 2000
        return out_val


class Output(InputOutput):
    def __init__(self, id):
        super().__init__(id)
        self.type = "OUTPUT"


class OutputOR(InputOutput):
    def __init__(self, gate_entry):
        id = gate_entry["identifier"]
        super().__init__(id)
        self.gate_entry = gate_entry
        self.type = "OUTPUT_OR2"
        pass

    def __call__(self, val1: float, val2: float, sim_settings) -> float:
        return val1 + val2


class OutputBuffer(InputOutput):
    def __init__(self, gate_entry):
        id = gate_entry["identifier"]
        super().__init__(id)
        self.gate_entry = gate_entry
        self.type = "OUTPUT_BUFFER"
        pass


class Gate(Device):
    def __init__(self, gate_entry):
        id = gate_entry["identifier"]
        super().__init__(id)

        self.gate_entry = gate_entry
        self.energy_rate = np.nan

        # self.identifier = gate_entry["identifier"]
        self.group = gate_entry["group"]

        pass

    @staticmethod
    def get_model_entry(gate_entry):
        model_entry = gate_entry["biorep"]["model"]
        return model_entry


class NOTGate(Gate):
    def __init__(self, gate_entry):
        super().__init__(gate_entry=gate_entry)
        self.type = "NOT"

        model_entry = Gate.get_model_entry(gate_entry)
        self.model = FourStatePromoterModel(model_entry)

        pass

    def __call__(self, in_val: float, sim_settings: dict) -> float:
        model = self.model

        mode = sim_settings["mode"]

        # ToDo Remove Conversion to RPU
        val = in_val  # * 2000

        mean_P, var_P, mean_M, var_M, energy_rate = model.get_distributions_and_energy_rate(val)

        # The energy is computed via the expected value
        self.energy_rate = energy_rate

        if mode == "samp":
            # The function value is sampled from the distribution
            common_val = var_P / mean_P ** 2 + 1
            mu = np.log(mean_P / np.sqrt(common_val))
            sigma = np.sqrt(np.log(common_val))

            sample = np.random.lognormal(mu, sigma, 1)

            output = sample[0]
        elif mode == "det":
            output = mean_P
        else:
            raise Exception(f"Mode {mode} is not supported.")

        # ToDo Remove Conversion to RPU
        out_val = output  # / 20000
        return out_val


class NORGate(Gate):

    def __init__(self, gate_entry):
        super().__init__(gate_entry=gate_entry)
        self.type = "NOR2"

        # self.or_gate = OutputOR(gate_entry)
        self.not_gate = NOTGate(gate_entry)
        pass

    def __call__(self, val1: float, val2: float, sim_settings: dict) -> float:
        # val = self.or_gate(val1, val2)
        val = val1 + val2
        output = self.not_gate(val, sim_settings=sim_settings)

        self.energy_rate = self.not_gate.energy_rate

        return output
