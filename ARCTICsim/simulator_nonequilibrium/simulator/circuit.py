from collections import OrderedDict

import numpy as np

from simulator.circuit_utils import CircuitAssignment, CircuitStructure
from simulator.libsim_coop import circuit_structure
from simulator.gatelib import GateLib
from simulator.utils import JsonFile


class Circuit:
    def __init__(self, structure: circuit_structure):
        self.structure = structure
        structure.combinational_order()

        self.propagation_graph = None
        self.inputs_graph = None
        self.generate_graph()

        self.assignment = None

        self.energy_rate = np.nan

        pass

    def generate_graph(self):

        structure = self.structure
        propagation_graph = OrderedDict()
        inputs_graph = OrderedDict()

        # {gate_id: iX for iX, gate_id in enumerate(self.structure.combinational_order())}

        for gate in structure.combinational_order():
            propagation_graph[gate] = [edge[1] for edge in structure.edges if edge[0] == gate]
            inputs_graph[gate] = [edge[0] for edge in structure.edges if edge[1] == gate]

        self.propagation_graph = propagation_graph
        self.inputs_graph = inputs_graph
        pass

    def set_assignment(self, assignment: CircuitAssignment):
        self.assignment = assignment
        self.energy_rate = np.nan
        pass




    def __call__(self, input_vals_dict):
        assignment = self.assignment
        propagation_graph = self.propagation_graph
        # inputs_graph = self.inputs_graph
        gate_input_vals = {id: [] for id in propagation_graph}
        gate_output_vals = {id: np.nan for id in propagation_graph}

        energy_rate = 0

        # Insert Input Values
        for input_id in input_vals_dict:
            gate_input_vals[input_id].append(input_vals_dict[input_id])

        # Propagate values through circuit
        for gate_id in propagation_graph:
            node_info = self.structure.node_infos[gate_id]
            device = assignment(node_info)
            str(device)
            in_vals = gate_input_vals[gate_id]
            assert not (any(np.isnan(in_vals)) or any([val < 0 for val in in_vals]))

            out_val = device(*in_vals)
            gate_output_vals[gate_id] = out_val

            # Propagate device output to subsequent gates
            for subsidary in propagation_graph[gate_id]:
                gate_input_vals[subsidary].append(out_val)

            #if node_info.type == "LOGIC":
            energy_rate += device.energy_rate

            #print("")
            #print(gate_id)
            #print(gate_input_vals)
            #print(gate_output_vals)
            #pass

        self.energy_rate = energy_rate

        return gate_output_vals


if __name__ == '__main__':
    json_structure = JsonFile(path="data/structures/structure_01110101.json")
    json_gatelib = JsonFile(path="../data/gate_libs/gate_lib_yeast.json")
    gatelib = GateLib(json_gatelib)

    # assignment_json = '{"a":"input_3","b":"input_1","c":"input_2","NOT_0":"H1_HlyIIR","NOT_2":"S2_SrpR","NOT_4":"B3_BM3R1","NOR2_1":"L1_LitR","NOR2_3":"P2_PhlF","O":"output_1"}'
    json_assignment = JsonFile(path="../data/assignments/assignment_01110101.json")
    assignment = CircuitAssignment(json_assignment.content, gate_lib=gatelib)

    structure = CircuitStructure(json_structure)

    circuit = Circuit(structure)
    circuit.set_assignment(assignment)
    print(circuit.propagation_graph)

    input_vals = {"input_1": 10 ** 9, "input_2": 10 ** 3, "input_3": 10 ** 8}
    input_vals = {"a": 0, "b": 1, "c": 0}
    output_vals = circuit(input_vals_dict=input_vals)
    print(output_vals)

    print(circuit.energy_rate)

    print(structure.to_dot())
    # N = 1000
    # out_vals = np.empty(N)
    # for iX in range(N):
    #     output_vals = circuit(input_vals_dict=input_vals)
    #     out_vals[iX] = output_vals["O"]


    pass
