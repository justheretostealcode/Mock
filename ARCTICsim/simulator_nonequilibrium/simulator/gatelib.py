import cProfile

from matplotlib import pyplot as plt

from simulator.particle_circuit_parts import ImplicitOrGate, NORGate, NOTGate, LutInput
from simulator.utils import JsonFile


# Class to represent the gatelib
class GateLib:
    def __init__(self, json_file: JsonFile):
        self.json = json_file

        if json_file is None:
            raise Exception("No Gate Lib information provided!")

        self.gates = []
        for gate_entry in json_file.data:
            gates = GateLib._populate_gate_entry(gate_entry)
            self.gates += gates

        self.gates_by_type_and_name = {}
        for gate in self.gates:
            prim_ident = gate.primitive_identifier
            ident = gate.identifier
            if prim_ident not in self.gates_by_type_and_name:
                self.gates_by_type_and_name[prim_ident] = {}

            if ident in self.gates_by_type_and_name[prim_ident]:
                raise Exception(f"{ident} is already in the lookup dict")

            self.gates_by_type_and_name[prim_ident][ident] = gate
        pass

    @staticmethod
    def _populate_gate_entry(gate_entry):
        def _get_gate_implementation(primitive_identifier):
            if primitive_identifier == "NOT":
                return NOTGate
            elif primitive_identifier == "NOR2":
                return NORGate
            elif primitive_identifier == "OR2":
                return ImplicitOrGate
            elif primitive_identifier == "INPUT":
                return LutInput
            else:
                raise Exception(f"Primitive Identifier \"{primitive_identifier}\" not supported")
            pass

        gates = []

        for primitive_identifier in gate_entry["primitiveIdentifier"]:
            gate_implementation = _get_gate_implementation(primitive_identifier)

            gate = gate_implementation(gate_entry)
            gates.append(gate)

        return gates



if __name__ == '__main__':
    gate_lib = GateLib(path="../data/gate_libs/gate_lib_yeast.json")

    not_1 = gate_lib.gates[1]

    profiler = cProfile.Profile()
    profiler.enable()

    for _ in range(10 ** 3):
        not_1(10 ** 8 * 1.0)

    profiler.disable()
    profiler.dump_stats("profile_results.prof")
    # profiler.print_stats()

    not_1 = gate_lib.gates[1]
    not_2 = gate_lib.gates[3]

    plt.figure()
    for not_i in [not_1, not_2]:
        samples = [not_i(10000000) for _ in range(1000)]
        plt.hist(samples, alpha=0.5, bins=100)

    plt.show()
    pass