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

    def __str__(self):
        return f"{self.__class__.__name__} ({self.identifier})"


class InputOutput(Device):
    def __call__(self, val):
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
        self.primitive_identifier = "INPUT"


        self.type = "INPUT"
        self.entry = input_entry

    def __call__(self, bool_val):
        return self.values[str(bool_val)]


class Output(InputOutput):
    def __init__(self, id):
        super().__init__(id)
        self.type = "OUTPUT"


class Gate(Device):
    def __init__(self, gate_entry):
        id = gate_entry["identifier"]
        super().__init__(id)

        self.gate_entry = gate_entry
        self.energy_rate = np.nan

        #self.identifier = gate_entry["identifier"]
        self.group = gate_entry["group"]

        pass

    @staticmethod
    def get_model_entry(gate_entry):
        model_entry = gate_entry["biorep"]["model"]
        return model_entry


class ImplicitOrGate(Gate):
    def __init__(self, gate_entry):
        super().__init__(gate_entry=gate_entry)
        self.primitive_identifier = "OR2"
        pass

    def __call__(self, val1, val2):
        self.energy_rate = 0
        return val1 + val2


class NOTGate(Gate):
    def __init__(self, gate_entry):
        super().__init__(gate_entry=gate_entry)
        self.primitive_identifier = "NOT"

        model_entry = Gate.get_model_entry(gate_entry)
        self.model = FourStatePromoterModel(model_entry)

        pass

    def __call__(self, val):
        model = self.model

        mean_P, var_P, mean_M, var_M, energy_rate = model.get_distributions_and_energy_rate(val)

        # The energy is computed via the expected value
        self.energy_rate = energy_rate

        # The function value is sampled from the distribution
        common_val = var_P / mean_P ** 2 + 1
        mu = np.log(mean_P / np.sqrt(common_val))
        sigma = np.sqrt(np.log(common_val))

        sample = np.random.lognormal(mu, sigma, 1)

        output = sample[0]
        return output


class NORGate(Gate):

    def __init__(self, gate_entry):
        super().__init__(gate_entry=gate_entry)
        self.primitive_identifier = "NOR2"

        self.or_gate = ImplicitOrGate(gate_entry)
        self.not_gate = NOTGate(gate_entry)
        pass

    def __call__(self, val1, val2):
        val = self.or_gate(val1, val2)
        output = self.not_gate(val)

        self.energy_rate = self.not_gate.energy_rate

        return output
