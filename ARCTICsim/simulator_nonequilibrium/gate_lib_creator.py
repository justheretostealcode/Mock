"""
Exemplary script to derive a Non-Equilibrium Compatible Library from the Cello for Yeast Library
"""
import numpy as np

from ARCTICsim.simulator_nonequilibrium.simulator.utils import JsonFile


def filter_by_collection(data_as_list, collection_name):
    filtered_data_as_list = filter(lambda elem: elem["collection"] == collection_name, data_as_list)
    return list(filtered_data_as_list)


def add_outputs(custom_lib, reporters):
    transcription_rate = 6
    translation_rate = 2
    rna_degradation_rate = 0.23104906018664842
    protein_degradation_rate = 0.0057762265
    cds_length = 561
    rna_length = 600
    protein_length = cds_length / 3

    for reporter in reporters:
        protein = {
            "identifier": "protein_" + reporter["name"],
            "name": reporter["name"],
            "collection": "proteins",
            "sequence_ids": [],
            "model_info": {
                "rna": {
                    "transcription_rate": transcription_rate,
                    "degradation_rate": rna_degradation_rate,
                    "e": 16,
                    "e_const": 0,
                    "length": rna_length
                },
                "protein": {
                    "translation_rate": translation_rate,
                    "degradation_rate": protein_degradation_rate,
                    "e": 42,
                    "e_const": 52,
                    "length": protein_length
                }
            }
        },
        device = {"identifier": "device_" + reporter["name"],
                  "name": reporter["name"],
                  "group": reporter["name"],
                  "collection": "devices",
                  "primitive_identifier": [
                      "OUTPUT_OR2",
                      "OUTPUT_BUFFER"
                  ],
                  "promoter_id": None,
                  "utr_id": None,
                  "cds_id": protein["identifier"],
                  "terminator_id": None
                  }

        custom_lib.append(protein)
        custom_lib.append(device)


def add_inputs(custom_lib, tf_inputs):
    # {"group": "aTc", "promoter": "ptet", "protein": "tetR", "low": 0.002, "high": 2.5},
    for tf_input in tf_inputs:
        protein = {"identifier": "tf_input_xylR",
                   "name": "xylR",
                   "collection": "tf_inputs",
                   "sequence_ids": [],
                   "model_info": {
                       "LUT": {
                           "0": 1,
                           "1": 0.01
                       }}}
        # ToDo match Promoter

        promoter = {"identifier": "promoter_" + tf_input["promoter"],
                    "name": tf_input["promoter"],
                    "collection": "promoters",
                    "cognate_transcription_factors": [
                        tf_input["protein"]
                    ],
                    "technology_mapping": None,
                    "model_info": {
                        "N_PROMOTER_STATES": 6,
                        "N_PARAMETERS": 14,
                        "INPUT_SCALING_FACTOR": 1,
                        "PROMOTER_ACTIVITY": [
                            0,
                            0,
                            0,
                            1,
                            1,
                            1
                        ],
                        "INFINITESIMAL_GENERATOR_FUNCTION": {
                            "expressions": [],
                            "matrices": {
                                "1": None,
                                "c": None,
                            }
                        }
                    },
                    "sequence_ids": []
                    }
        device = {"identifier": "input_" + tf_input["group"],
                  "name": tf_input["group"],
                  "group": tf_input["group"],
                  "collection": "devices",
                  "primitive_identifier": [
                      "INPUT"
                  ],
                  "promoter_id": promoter["identifier"],
                  "utr_id": None,
                  "cds_id": protein["identifier"],
                  "terminator_id": None
                  }


        custom_lib.append(protein)
        custom_lib.append(promoter)
        custom_lib.append(device)
        pass


def match_promoter(response_characteristic:dict):
    """
    This method performs the parameter matching for the promoter model to the response characteristics

    :param response_characteristic: dict providing the input values and corresponding outputcharacteristics
    :return:
    """

    # Steps
    # Setup Model to train (Gene includes YFP output and the input is a TF Level (possibly derived from RPU))
    # Match model to training data



if __name__ == '__main__':
    # The cello for yeast library can be downloaded in
    # the Supplementary Information section of https://www.nature.com/articles/s41564-020-0757-2#additional-information
    cello_library_path = "ARCTICsim/simulator_nonequilibrium/data/reference_data/yeast/SC1C1G1T1.UCF.json"
    cello_library_json = JsonFile(path=cello_library_path)
    cello_library = cello_library_json.data
    collections = set([elem["collection"] for elem in cello_library])

    cello_models = filter_by_collection(cello_library, "models")
    cello_functions = filter_by_collection(cello_library, "functions")
    cello_structures = filter_by_collection(cello_library, "structures")
    cello_cytometry = list(filter(lambda elem: elem["name"].endswith("_cytometry"), cello_functions))
    cello_toxicity = list(filter(lambda elem: elem["name"].endswith("_toxicity"), cello_functions))
    cello_parts = filter_by_collection(cello_library, "parts")

    cello_library_by_collection = {collection: [elem for elem in cello_library if elem["collection"] == collection]
                                   for
                                   collection in collections}
    cello_structures_by_name = {elem["name"]: elem for elem in cello_structures}
    cello_parts_by_name = {elem["name"]: elem for elem in cello_parts}
    cello_cytometry_by_name = {elem["name"]: elem for elem in cello_cytometry}
    cello_toxicity_by_name = {elem["name"]: elem for elem in cello_cytometry}

    # Useful for conversion to other gatelibrary formats such as the cello E. Coli lib.
    dict_ids = {"SPECIES": "Yeast",
                "name": "name",
                "table": "table",
                "bin": "bin",
                "output": "output",
                "x": "x"}

    """
    The relevant collections are
    "gates":        Includes an overview on all the genetic gates present, their names,
                    associated structures, transcription factors and gate types.
    "models":       Includes the cello model representation of the genetic gates
    "structures":   Defines the the genetic gates with output devices (promoters) 
                    and corresponding inner signaling molecules (transcription factors). 
                    Also includes the codon variants of the promoters as well as the 
                    proteins and the kozak sequences to use.
    "functions":    The information included are input to output relations. Of particular 
                    importance is the cytometry data identifiable with the suffix 
                    "_cytometry" preced by the associated genetic gate's name. Also the 
                    characteristics of the toxicity are captured (suffix "_toxicity") and
                    equations for the "Hill_response" and "linear_input_composition" 
                    provided. 
    "parts":        Provides the DNA sequence of all components. 
                    
    Conventions     
                    Naming format of genetic gates in the cello library
                    The genetic gates are indicated by a name including an underscore, for 
                    example "P1_LexA". The part following the underscore gives rise to the 
                    transcription factor used inside the gate. With this information, the 
                    first part then determines the promoter used. In the context of the example
                    "P1" refers to "pLexA1" as it is the first variant of the promoter for the 
                    transcription factor "LexA". The second promoter variant is "pLexA2"
                    Further underscores add information on the codon variant (e.g. "pLexA1_a" or 
                    "pLexA1_b") or characterize the information provided for the genetic gate as
                    "P1_LexA_structure" identifies the structure of the genetic gate "P1_LexA",
                    while "P1_LexA_toxicity" gives rise to the corresponding toxicity information.
                    
                    
                    
    
    """

    cello_gates = cello_library_by_collection["gates"]
    custom_lib = []

    # Reporter Output
    reporters = [{"name": "YFP"}]
    add_outputs(custom_lib, reporters=reporters)


    # Transcription Factor Inputs
    tf_inputs = [{"group": "aTc", "promoter": "ptet", "protein": "tetR", "low": 0.002, "high": 2.5},
                 {"group": "xyl", "promoter": "pxyl", "protein": "xylR", "low": 0.003, "high": 1.8},
                 {"group": "IPTG", "promoter": "plac", "protein": "lacl", "low": 0.01, "high": 1.3}]

    add_inputs(custom_lib, tf_inputs=tf_inputs)


    # Actual gates
    for cello_gate in cello_gates:
        gate_name = cello_gate["name"]
        structure = cello_structures_by_name[cello_gate["structure"]]

        promoter_name = structure["outputs"][0].split("_")[0]
        promoter_sequence_names = structure["outputs"]
        utr_name = [elem for elem in structure["devices"] if "cassette" in elem["name"]][0]["components"][0]
        protein_name = cello_gate["regulator"]
        protein_sequence_names = [elem["components"][1] for elem in structure["devices"] if
                                  "cassette" in elem["name"]]

        cello_utr_part = cello_parts_by_name[utr_name]
        cello_promoter_part = [cello_parts_by_name[name] for name in promoter_sequence_names]
        cello_regulator_parts = [cello_parts_by_name[name] for name in protein_sequence_names]

        # protein_sequences = [elem for elem in structure["devices"] if "cassette" in elem["name"]][0]["components"][0]

        # associated_parts = [elem for elem in cello_parts if elem["name"] in ]

        # ToDo derive parameters of promoter

        utr_sequence = {"collection": "sequences",
                        "identifier": "sequence_" + utr_name,
                        "name": utr_name,
                        "sequence": cello_utr_part["dnasequence"]}
        protein_sequences = [{"collection": "sequences",
                              "identifier": "sequence_" + name,
                              "name": name,
                              "sequence": cello_parts_by_name[name]["dnasequence"]}
                             for name in protein_sequence_names]
        promoter_sequences = [{"collection": "sequences",
                               "identifier": "sequence_" + name,
                               "name": name,
                               "sequence": cello_parts_by_name[name]["dnasequence"]}
                              for name in promoter_sequence_names]

        utr = {"collection": "utrs",
               "identifier": "utr_" + utr_name,
               "name": utr_name,
               "sequence_id": utr_sequence["identifier"]}

        cds_length = np.median([len(seq["sequence"]) for seq in protein_sequences])
        rna_length = len(utr_sequence["sequence"]) + cds_length
        protein_length = cds_length / 3

        # ToDo Redefine rates
        transcription_rate = 6
        translation_rate = 2
        rna_degradation_rate = 0.23104906018664842
        protein_degradation_rate = 0.0057762265
        protein = {"identifier": "protein_" + protein_name,
                   "name": protein_name,
                   "collection": "proteins",
                   "sequence_ids": [seq["identifier"] for seq in protein_sequences],
                   "model_info": {
                       "rna": {
                           "transcription_rate": transcription_rate,
                           "degradation_rate": rna_degradation_rate,  # ToDo Redefine rates
                           "e": 16,
                           "e_const": 0,
                           "length": rna_length,
                       },
                       "protein": {
                           "translation_rate": translation_rate,  # ToDo Redefine rates
                           "degradation_rate": protein_degradation_rate,  # ToDo Redefine rates
                           "e": 42,
                           "e_const": 52,
                           "length": protein_length,
                       }
                   }}
        cognate_transcription_factors = [protein_name, protein_sequence_names]
        promoter = {"identifier": "promoter_" + promoter_name,
                    "name": promoter_name,
                    "collection": "promoters",
                    "cognate_transcription_factors": cognate_transcription_factors,
                    "sequence_ids": [seq["identifier"] for seq in promoter_sequence_names],
                    "technology_mapping": None,
                    "model_info": {
                        "N_PROMOTER_STATES": 6,
                        "N_PARAMETERS": 14,
                        "INPUT_SCALING_FACTOR": 0.0001,
                        "PROMOTER_ACTIVITY": [0, 0, 0, 1, 1, 1],
                        "INFINITESIMAL_GENERATOR_FUNCTION": {"matrices": {"1": [[]],
                                                                          "c": [[]]
                                                                          }}},
                    }

        device = {"collection": "devices",
                  "identifier": "device_" + gate_name,
                  "name": gate_name,
                  "group": cello_gate["group"],
                  "regulator": cello_gate["regulator"],
                  "primitive_identifier": ["NOT",
                                           "NOR2"],
                  "utr_id": utr["identifier"],  # ToDo Add the corresponding ID
                  "cds_id": protein["identifier"],  # ToDo Add the corresponding ID
                  "terminator_id": None,  # terminators are currently not included in the library
                  "promoter_id": promoter["identifier"],
                  "color": cello_gate["color"]
                  }  # ToDo Add the corresponding ID

        # Next steps:
        # Match Characteristics of promoter to cytometry data
        # Infer Inflection Point
        # Store everything in a list

        # Add everything to the lists

        pass
        pass
        # Write out everything to a json file
