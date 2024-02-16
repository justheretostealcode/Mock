import itertools
import json
import numpy as np
from gvgen import GvGen

from ARCTICsim.simulator_nonequilibrium.simulator.gatelib import GateLib, Device
from ARCTICsim.simulator_nonequilibrium.simulator.utils import JsonFile
from ARCTICsim.simulator_thermodynamics.libsim import circuit_structure


#####################################################################################################
#                                                                                                   #
# The subsequent classes represent extensions to classes in libsim_coop.py from Nicolai Engelmann   #
#                                                                                                   #
#####################################################################################################

class CircuitStructure(circuit_structure):
    class NodeInfo:
        def __init__(self, node_info_dict):
            self.type = node_info_dict["type"]
            # Legacy fields
            # self.primitiveIdentifier = node_info_dict["primitiveIdentifier"]
            # self.expression = node_info_dict["expression"]
            self.id = node_info_dict["id"]

    def __init__(self, json: JsonFile):
        super().__init__(json.content)

        self.node_infos = {node_info["id"]: CircuitStructure.NodeInfo(node_info) for node_info in
                           json.data['graph']['nodes']}

        self.truthtable = TruthTable(json.data["truthtable"])
        pass

    def get_outgoing_edges(self, node_id):
        return [edge for edge in self.edges if edge[0] == node_id]

    def get_ingoing_edges(self, node_id):
        return [edge for edge in self.edges if edge[1] == node_id]

    def to_dot(self):
        g = GvGen()

        items = {}
        for node in self.nodes:
            item = g.newItem(node)
            items[node] = item

        for edge in self.edges:
            src_item = items[edge[0]]
            dest_item = items[edge[1]]
            g.newLink(src_item, dest_item)

        dot_rep = g.dot()
        return dot_rep


#####################################################################################################
#                                                                                                   #
# Own Classes                                                                                       #
#                                                                                                   #
#####################################################################################################


class CircuitAssignment:
    def __init__(self, json_file: JsonFile, gate_lib: GateLib):
        self.mapping = json_file.data
        self.gate_lib = gate_lib

    def __call__(self, node_info: CircuitStructure.NodeInfo) -> Device:

        id = node_info.id
        # Legacy field
        # primitiveIdentifier = node_info.primitiveIdentifier
        node_type = node_info.type

        device_id = self.mapping[id]

        raise Exception("Not Compatible Yet: Adapt to new GateLib")
        # ToDo Adapt to new GateLib
        gates_by_type_and_name = self.gate_lib.gates_by_type_and_name
        if node_type not in gates_by_type_and_name:
            raise Exception(
                f"The node type {node_type} is not supported. It needs to be either of {self.gate_lib.gates_by_type_and_name.keys()}.")

        if device_id not in gates_by_type_and_name[node_type]:
            raise Exception(
                f"The device {device_id} is not known. The known devices for node type {node_type} are {self.gate_lib.gates_by_type_and_name[node_type].keys()}.")

        gate = self.gate_lib.gates_by_type_and_name[node_type][device_id]
        return gate


class TruthTable:
    def __init__(self, truthtable_string):
        circ_outs_reversed = list(map(int, truthtable_string))
        circ_outs = circ_outs_reversed[::-1]

        num_inputs = int(np.log2(len(circ_outs)))
        assert 2 ** num_inputs == len(circ_outs)

        truthtable = np.zeros(shape=(len(circ_outs), num_inputs + 1), dtype=int)

        truthtable[:, 0:3] = np.array(list(itertools.product([0, 1], [0, 1], [0, 1])))
        truthtable[:, 3] = circ_outs

        self.num_inputs = num_inputs
        self.truthtable = truthtable

    def input_output_truthtable(self):
        input_output_truthtable = [(row[0:3], row[3]) for row in self.truthtable]
        return input_output_truthtable
        # return self.truthtable[:, 0:3], self.truthtable[:, 3]

    def __str__(self):
        string_rep = ""
        truthtable = self.input_output_truthtable()

        for entry in truthtable:
            string_rep += str(entry[0]) + "|" + str(entry[1]) + "\n"

        return string_rep


def load_structure(settings: dict) -> CircuitStructure:
    structure_path = settings["structure"]
    json_structure = JsonFile(path=structure_path)
    structure = CircuitStructure(json_structure)
    return structure
