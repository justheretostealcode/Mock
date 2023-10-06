import json

from ARCTICsim.adaptation.data_structures import ChromosomalLocation, Plasmid, Position
from ARCTICsim.simulator_nonequilibrium.simulator.circuit_utils import CircuitStructure, CircuitAssignment
from ARCTICsim.simulator_nonequilibrium.simulator.gatelib import GateLib
from ARCTICsim.simulator_nonequilibrium.simulator.particle_circuit_parts import LutInput, InputOutput, Output
from ARCTICsim.simulator_nonequilibrium.simulator.utils import JsonFile


coding_sequences = {}
def get_coding_sequence(transcription_factor):
    if transcription_factor not in coding_sequences:
        coding_sequences[transcription_factor] = {"FREE": [transcription_factor + "_" + suffix for suffix in ["b", "a"]], "USED": []}

    cds = coding_sequences[transcription_factor]["FREE"].pop()
    coding_sequences[transcription_factor]["USED"].append(cds)

    return cds

promoters = {}
def get_promoter(promoter):
    if promoter not in promoters:
        promoters[promoter] = {"FREE": [promoter + "_" + suffix for suffix in ["b", "a"]], "USED": []}

    promoter_name = promoters[promoter]["FREE"].pop()
    promoters[promoter]["USED"].append(promoter_name)

    return promoter_name

def place_on_plasmid(structure, assignment):
    def get_transcription_factor_and_promoter(genetic_gate):
        tf = genetic_gate.transcription_factor
        promoter = genetic_gate.promoter

        if tf is None or promoter is None:
            id = genetic_gate.identifier

            promoter2 = id

            if isinstance(genetic_gate, Output):
                tf2 = "YFP"
            elif isinstance(genetic_gate, InputOutput):
                # ToDo Add Group
                tf2 = None
            else:
                promoter2, tf2 = id.split("_")
                promoter2 = id

            if tf is None:
                tf = tf2

            if promoter is None:
                promoter = promoter2

        return tf, promoter

    plasmid = Plasmid()
    # chromosomal_location = ChromosomalLocation()
    nodes = structure.combinational_order()

    plasmid_pos = 0

    for node_id in nodes:
        node_info = structure.node_infos[node_id]
        source_gate = assignment(node_info=node_info)
        source_tf, source_promoter = get_transcription_factor_and_promoter(source_gate)
        source_promoter = get_promoter(source_promoter)

        outgoing_edges = structure.get_outgoing_edges(node_id=node_id)
        target_nodes = [edge[1] for edge in outgoing_edges]
        target_gates = [assignment(node_info=structure.node_infos[target_node]) for target_node in target_nodes]
        for target_gate in target_gates:
            target_tf, target_promoter = get_transcription_factor_and_promoter(target_gate)
            chromosomal_id = plasmid_pos % 2
            chromosomal_location = plasmid_pos // 2

            target_cds = get_coding_sequence(target_tf)
            position = Position(id=chromosomal_location, promoter=source_promoter, coding_sequence=target_cds)
            plasmid[chromosomal_id][chromosomal_location] = position

            #location.promoter = source_promoter
            #location.coding_sequence = target_tf

            plasmid_pos += 1
        pass

    return plasmid

def to_json(data, file_path):
    with open(file_path, "w") as file:
        json.dump(data, file, indent=4)

def from_json(file_path):
    with open(file_path, "r") as file:
        data = json.load(file)
    return data

if __name__ == '__main__':
    structure_path = "ARCTICsim/simulator_nonequilibrium/data/structures/structure_01110101.json"
    assignment_str = '{"a":"input_3","b":"input_1","c":"input_2","NOT_0":"P1_PsrA","NOT_2":"P1_IcaR","NOT_4":"P1_PhlF","NOR2_1":"P1_QacR","NOR2_3":"P1_HKCI","O":"output_1"}'
    gatelib_path = "ARCTICsim/simulator_nonequilibrium/data/gate_libs/gate_lib_yeast_cello_data_fixed_limits.json"

    json_structure = JsonFile(path=structure_path)
    json_gatelib = JsonFile(path=gatelib_path)
    json_assignment = JsonFile(content=assignment_str)

    structure = CircuitStructure(json=json_structure)
    gate_lib = GateLib(json_file=json_gatelib)
    assignment = CircuitAssignment(json_file=json_assignment, gate_lib=gate_lib)

    plasmid = place_on_plasmid(structure, assignment)

    print(assignment.mapping)

    print(plasmid)

    plasmid_json = plasmid.to_json()
    to_json(plasmid_json, "plasmid.json")
    # plasmid_json_2 = from_json("plasmid.json")
    # plasmid_2 = Plasmid.from_json(plasmid_json_2)
    # plasmid_json_2 = plasmid_2.to_json()
    # to_json(plasmid_json_2, "plasmid_2.json")


    raise Exception("Adapt to the use of multiple promoters, if promoter is used multiple times. Identify the meaning of P1 und P2 for this.")
    raise Exception("Output as JSON")
    pass
