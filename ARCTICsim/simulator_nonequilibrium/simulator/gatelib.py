import cProfile

import numpy as np
from matplotlib import pyplot as plt

from ARCTICsim.simulator_nonequilibrium.models.moment_model_new import RNAMomentModel, ProteinMomentModel, \
    CombinedMomentModel
from ARCTICsim.simulator_nonequilibrium.models.promoter_model_new import PromoterModel
from ARCTICsim.simulator_nonequilibrium.simulator.particle_circuit_parts import NOTGate, NORGate, LutInput, OutputOR, \
    OutputBuffer
from ARCTICsim.simulator_nonequilibrium.simulator.utils import JsonFile


# Class to represent the gatelib
class GateLib:
    def __init__(self, json_file: JsonFile):
        self.json = json_file

        if json_file is None:
            raise Exception("No Gate Library information provided!")

        self.gates = []
        for gate_entry in json_file.data:
            gates = GateLib._populate_gate_entry(gate_entry)
            self.gates += gates

        self.gates_by_type_and_name = {}
        for gate in self.gates:
            gate_type = gate.type
            ident = gate.identifier
            if gate_type not in self.gates_by_type_and_name:
                self.gates_by_type_and_name[gate_type] = {}

            if ident in self.gates_by_type_and_name[gate_type]:
                raise Exception(f"{ident} is already in the lookup dict")

            self.gates_by_type_and_name[gate_type][ident] = gate
        pass

    @staticmethod
    def _populate_gate_entry(gate_entry):
        def _get_gate_implementation(gate_type):
            if gate_type == "NOT":
                return NOTGate
            elif gate_type == "NOR2":
                return NORGate
            elif gate_type == "INPUT":
                return LutInput
            elif gate_type == "OUTPUT_OR2":
                return OutputOR
            elif gate_type == "OUTPUT_BUFFER":
                return OutputBuffer
            else:
                raise Exception(f"Type \"{gate_type}\" not supported")
            pass

        gates = []

        for gate_type in gate_entry["primitiveIdentifier"]:
            gate_implementation = _get_gate_implementation(gate_type)

            gate = gate_implementation(gate_entry)
            gates.append(gate)

        return gates


class GateLibCollectionBased:
    def __init__(self, json_file: JsonFile):
        self.json = json_file

        if json_file is None:
            raise Exception("No Gate Library information provided!")

        """
        Gate Lib Format
        The new gate lib design is collection based
        
        The collections included are:
        - Promoter (and associated TF)
        - Transcription Factors
        - Sequences (Promoter and Transcription Factor)
        - Devices (Pairs of TF and Promoter with associated Logic Type Realized, Includes OutputBuffer and OutputOR2, LUT Inputs)
        
        
        Comment:
        As the inputs represent promoter activity, the LUT Inputs fall into the same category as the promoters.
        """

        # self.gates = []
        # for gate_entry in json_file.data:
        #     gates = GateLib._populate_gate_entry(gate_entry)
        #     self.gates += gates
        #
        # self.gates_by_type_and_name = {}
        # for gate in self.gates:
        #     gate_type = gate.type
        #     ident = gate.identifier
        #     if gate_type not in self.gates_by_type_and_name:
        #         self.gates_by_type_and_name[gate_type] = {}
        #
        #     if ident in self.gates_by_type_and_name[gate_type]:
        #         raise Exception(f"{ident} is already in the lookup dict")
        #
        #     self.gates_by_type_and_name[gate_type][ident] = gate
        pass


class GateLibEntry:
    def __init__(self, collection_identifier: str = None, gate_lib_entry_dict: dict = None):
        if collection_identifier is None:
            raise Exception("Collection Identifier is not set!")

        if gate_lib_entry_dict is None:
            raise Exception("No gate_lib_entry_dict provided.")

        self.collection_identifier = collection_identifier
        self.gate_lib_entry_dict = gate_lib_entry_dict
        self.id = gate_lib_entry_dict["identifier"]
        self.collection = gate_lib_entry_dict["collection"]
        self.group = gate_lib_entry_dict["group"] if "group" in gate_lib_entry_dict else None

        if self.collection != self.collection_identifier:
            raise Exception(
                f"Mismatch in collection (Trying to instantiate {self.collection_identifier} with an entry of collection {self.collection}")

    def __str__(self):
        return f"{type(self)} {self.id}"


class GateLibModelEntry(GateLibEntry):
    def __init__(self, collection_identifier: str = None, gate_lib_entry_dict: dict = None):
        super().__init__(collection_identifier=collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)

        self.model_info = gate_lib_entry_dict["model"]
        self.model = self._populate_model()

    def _populate_model(self):
        """
        Is implemented by the respective subclass to realize the model creation.
        :return: The model populated with the parameters provided in self.model_info
        """
        raise Exception("Not Implemented Yet")

    def __call__(self, in_vals: dict, sim_settings: dict, *args, **kwargs):
        """
        Is implemented by the respective subclass or the subclass provides a model matching the interface.
        :param args:
        :param kwargs:
        :return:
        """
        output = self.model(in_vals, sim_settings=sim_settings)
        return output


class Promoter(GateLibModelEntry):
    def __init__(self, gate_lib_entry_dict: dict):
        self.collection_identifier = "promoters"
        super().__init__(collection_identifier=self.collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)

        self.cognate_transcription_factors = gate_lib_entry_dict["cognate_transcription_factors"]
        self.sequence_ids = gate_lib_entry_dict["sequence_ids"] # A single construct can be asssociated to multiple sequences (for example codon alternatives)

    def _populate_model(self):
        model_info = self.model_info
        num_states = model_info["N_PROMOTER_STATES"]
        num_trainable_parameters = model_info["N_PARAMETERS"]
        promoter_activity = model_info["PROMOTER_ACTIVITY"]
        infinitesimal_generator_function = model_info["INFINITESIMAL_GENERATOR_FUNCTION"]

        model = PromoterModel(num_states,
                              num_trainable_parameters,
                              infinitesimal_generator_function,
                              per_state_promoter_activity=promoter_activity)

        return model

    # def __call__(self, cell_state: dict, sim_settings: dict, *args, **kwargs):
    #     # Single input promoter
    #     # c subsumes cognate TFs
    #
    #     output = super()(in_vals, sim_settings)
    #     return output


class LutPromoter(GateLibModelEntry):
    def __init__(self, gate_lib_entry_dict: dict):
        self.collection_identifier = "lut_promoters"
        super().__init__(collection_identifier=self.collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)

        self.cognate_transcription_factors = gate_lib_entry_dict["cognate_transcription_factors"]
        # The corresponding sequences are currently not included in the library.
        self.sequence_idds = []

    def _populate_model(self):
        model_info = self.model_info
        LUT = {float(key): model_info["LUT"][key] for key in model_info["LUT"]}
        # ToDo Check whether this works as intended
        def model(in_val, sim_settings=None): LUT[in_val]

        return model


class CodingSequence(GateLibModelEntry):
    """
    Superclass of expresseble sequences.
    """

    def __init__(self, gate_lib_entry_dict: dict):
        self.collection_identifier = "coding_sequence"
        super().__init__(collection_identifier=self.collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)

        self.sequence_ids = gate_lib_entry_dict["sequence_ids"] # A single construct can be asssociated to multiple sequences (for example codon alternatives)



class Protein(CodingSequence):
    def __init__(self, gate_lib_entry_dict: dict):
        self.collection_identifier = "protein"
        super().__init__(collection_identifier=self.collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)

        self.name = gate_lib_entry_dict["name"]
    def _populate_model(self):
        model_info = self.model_info
        # Includes molecule dependent transcription, RNA degradation and the respective energy requirements and lengths.

        transcription_rate = model_info["transcription_rate"]
        rna_degradation_rate = model_info["degradation_rate"]
        energy_per_nucleotide = model_info["energy_per_nucleotide"]
        energy_per_rna = model_info["energy_per_rna"]
        rna_length = model_info["length"]

        translation_rate = model_info["translation_rate"]
        protein_degradation_rate = model_info["degradation_rate"]
        energy_per_amino_acid = model_info["energy_per_amino_acid"]
        energy_per_protein = model_info["energy_per_protein"]
        protein_length = model_info["length"]

        model = CombinedMomentModel(transcription_rate, rna_degradation_rate, energy_per_nucleotide, energy_per_rna,
                                    rna_length,
                                    translation_rate, protein_degradation_rate, energy_per_amino_acid,
                                    energy_per_protein, protein_length)

        return model


# class RNADynamics(GateLibModelEntry):
#     def __init__(self, gate_lib_entry_dict):
#         self.collection_identifier = "rna_dynamics"
#         super().__init__(collection_identifier=self.collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)
#
#     def _populate_model(self):
#         model_info = self.model_info
#         # Includes molecule dependent transcription, RNA degradation and the respective energy requirements and lengths.
#
#         transcription_rate = model_info["transcription_rate"]
#         degradation_rate = model_info["degradation_rate"]
#         energy_per_nucleotide = model_info["energy_per_nucleotide"]
#         energy_per_rna = model_info["energy_per_rna"]
#         length = model_info["length"]
#         model = RNAMomentModel(transcription_rate, degradation_rate,
#                                energy_per_nucleotide, energy_per_rna,
#                                length)
#         return model
#
#
# class ProteinDynamics(GateLibModelEntry):
#     def __init__(self, gate_lib_entry_dict):
#         self.collection_identifier = "protein_dynamics"
#         super().__init__(collection_identifier=self.collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)
#
#     def _populate_model(self):
#         model_info = self.model_info
#         # Includes molecule dependent translation, protein degradation and the respective energy requirements and lengths.
#
#         translation_rate = model_info["translation_rate"]
#         degradation_rate = model_info["degradation_rate"]
#         energy_per_amino_acid = model_info["energy_per_amino_acid"]
#         energy_per_protein = model_info["energy_per_protein"]
#         length = model_info["length"]
#         model = ProteinMomentModel(translation_rate, degradation_rate,
#                                    energy_per_amino_acid, energy_per_protein,
#                                    length)
#         return model


# class Protein(GateLibModelEntry):
#     def __init__(self, gate_lib_entry_dict):
#         self.collection_identifier = "proteins"
#         super().__init__(collection_identifier=self.collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)
#
#     def _populate_model(self):
#         raise Exception("Not Implemented Yet")


class Sequence(GateLibEntry):
    def __init__(self, gate_lib_entry_dict: dict):
        self.collection_identifier = "sequences"
        super().__init__(collection_identifier=self.collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)

        self.type = gate_lib_entry_dict["type"]
        self.sequence = gate_lib_entry_dict["dna_sequence"]


class UTR(GateLibEntry):
    def __init__(self, gate_lib_entry_dict: dict):
        self.collection_identifier = "utrs"
        super().__init__(collection_identifier=self.collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)


class Device(GateLibEntry):
    """
    Represents the "Genetic Gate" in the meaning of Nielsen et al. 2016

    Needs to be populated last.
    """

    def __init__(self, gate_lib_entry_dict: dict, parts: dict):
        self.collection_identifier = "devices"
        super().__init__(collection_identifier=self.collection_identifier, gate_lib_entry_dict=gate_lib_entry_dict)

        self.promoter_id = gate_lib_entry_dict["promoter_id"]
        self.utr_id = gate_lib_entry_dict["utr_id"]
        self.cds_id = gate_lib_entry_dict["cds_id"]
        self.terminator_id = gate_lib_entry_dict["terminator_id"]
        self.ids = [self.promoter_id, self.utr_id, self.cds_id, self.terminator_id]
        if parts is not None:
            parts_to_use = dict(parts)
            parts_to_use = {key: parts_to_use[key] if key in parts_to_use else None for key in self.ids}

            self.promoter = parts_to_use[self.promoter_id]
            self.utr = parts_to_use[self.utr_id]
            self.cds = parts_to_use[self.cds_id]
            self.terminator_id = parts_to_use[self.terminator_id]

    def __str__(self):
        return "Device: " + ", ".join(self.ids)


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
