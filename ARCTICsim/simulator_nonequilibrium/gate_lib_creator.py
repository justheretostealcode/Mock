"""
Exemplary script to derive a Non-Equilibrium Compatible Library from the Cello for Yeast Library
"""
import json

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize

from ARCTICsim.simulator_nonequilibrium.simulator.circuit import Gene
from ARCTICsim.simulator_nonequilibrium.simulator.gatelib import Promoter, UTR, Protein, TranscriptionfactorInputs, \
    SensorPromoter
from ARCTICsim.simulator_nonequilibrium.simulator.histogram_data import GateCytometryData
from ARCTICsim.simulator_nonequilibrium.simulator.utils import JsonFile


def filter_by_collection(data_as_list, collection_name):
    filtered_data_as_list = filter(lambda elem: elem["collection"] == collection_name, data_as_list)
    return list(filtered_data_as_list)


def get_promoter_entry(name, cognate_transcription_factors=None, sequence_ids=None, scaling_factor=1):
    if cognate_transcription_factors is None:
        cognate_transcription_factors = []
    if sequence_ids is None:
        sequence_ids = []
    promoter_entry = {"identifier": "promoter_" + name,
                      "name": name,
                      "collection": "promoters",
                      "cognate_transcription_factors": cognate_transcription_factors,
                      "sequence_ids": sequence_ids,
                      "technology_mapping": None,
                      "model_info": {"N_PROMOTER_STATES": 6,
                                     "N_PARAMETERS": 16,
                                     "INPUT_SCALING_FACTOR": scaling_factor,
                                     "PROMOTER_ACTIVITY": [0, 0, 0, 1, 1, 1],
                                     "INFINITESIMAL_GENERATOR_FUNCTION": {"matrices": {"1": [[]],
                                                                                       "c": [[]]
                                                                                       }}}
                      }

    return promoter_entry


def get_protein_entry(name, sequence_ids=None,
                      transcription_rate=None, rna_degradation_rate=None, rna_length=None,
                      translation_rate=None, protein_degradation_rate=None, protein_length=None):
    if sequence_ids is None:
        sequence_ids = []
    protein_entry = {"identifier": "protein_" + name,
                     "name": name,
                     "collection": "proteins",
                     "sequence_ids": sequence_ids,
                     "model_info": {"rna": {"transcription_rate": transcription_rate,
                                            "degradation_rate": rna_degradation_rate,
                                            "e": 16,
                                            "e_const": 0,
                                            "length": rna_length
                                            },
                                    "protein": {"translation_rate": translation_rate,
                                                "degradation_rate": protein_degradation_rate,
                                                "e": 42,
                                                "e_const": 52,
                                                "length": protein_length
                                                }
                                    }
                     }
    return protein_entry


def get_device_entry(name, group, regulator=None, primitive_identifer=None, utr_id=None,
                     cds_id=None, promoter_id=None, terminator_id=None, color=None):
    if primitive_identifer is None:
        primitive_identifer = []

    device_entry = {"collection": "devices",
                    "identifier": "device_" + name,
                    "name": name,
                    "group": group,
                    "regulator": regulator,
                    "primitive_identifier": primitive_identifer,
                    "utr_id": utr_id,  # ToDo Add the corresponding ID
                    "cds_id": cds_id,  # ToDo Add the corresponding ID
                    "terminator_id": terminator_id,  # terminators are currently not included in the library
                    "promoter_id": promoter_id,
                    "color": color
                    }  # ToDo Add the corresponding ID
    return device_entry


def get_sequence_entry(name, sequence=""):
    sequence_entry = {"collection": "sequences",
                      "identifier": "sequence_" + name,
                      "name": name,
                      "sequence": sequence}
    return sequence_entry


def add_outputs(custom_lib, reporters, rates:dict):
    transcription_rate = rates["transcription_rate"]
    translation_rate = rates["translation_rate"]
    rna_degradation_rate = rates["rna_degradation_rate"]
    protein_degradation_rate = rates["protein_degradation_rate"]
    cds_length = 561
    rna_length = 600
    protein_length = cds_length / 3
    for reporter in reporters:
        utr_entry = {"collection": "utrs",
                     "identifier": "utr_" + reporter["utr_name"],
                     "name": reporter["utr_name"],
                     "sequence_id": None}

        protein_entry = get_protein_entry(name=reporter["name"],
                                          sequence_ids=[],
                                          transcription_rate=transcription_rate,
                                          rna_degradation_rate=rna_degradation_rate,
                                          rna_length=rna_length,
                                          translation_rate=translation_rate,
                                          protein_degradation_rate=protein_degradation_rate,
                                          protein_length=protein_length)

        device_entry = get_device_entry(name=reporter["name"], group=reporter["name"],
                                        primitive_identifer=["OUTPUT_OR2", "OUTPUT_BUFFER"],
                                        utr_id=utr_entry["identifier"],
                                        cds_id=protein_entry["identifier"])

        custom_lib.append(utr_entry)
        custom_lib.append(protein_entry)
        custom_lib.append(device_entry)


def add_inputs(custom_lib, tf_inputs, reporter_information):
    # {"group": "aTc", "promoter": "ptet", "protein": "tetR", "low": 0.002, "high": 2.5},
    reference_sensor_promoter_entry = None
    for tf_input in tf_inputs:
        # protein_entry = {"identifier": "tf_input_" + tf_input["protein"],
        #                  "name": tf_input["protein"],
        #                  "collection": "tf_inputs",
        #                  "sequence_ids": [],
        #                  "model_info": {
        #                      "LUT": {
        #                          "0": 10,
        #                          "1": 0.001
        #                      }}}
        # ToDo match Promoter
        promoter_entry = get_promoter_entry(name=tf_input["promoter"],
                                            cognate_transcription_factors=[tf_input["group"], tf_input["protein"]],
                                            sequence_ids=None)
        promoter_entry["collection"] = "sensor_promoters"
        promoter_entry["identifier"] = "sensor_" + promoter_entry["identifier"]
        data = tf_input["data"]
        promoter_entry["model_info"] = {"low": tf_input["low"],
                                        "high": tf_input["high"],
                                        "LUT": data}
        promoter_entry["technology_mapping"] = {"y_min": data[tf_input["low"]],
                                                "y_max": data[tf_input["high"]]}

        # response_characteristic = {"X": np.array(list(data.keys())),
        #                            "Y": np.array([[data[inducer_concentration],
        #                                            0.25 * data[inducer_concentration]
        #                                            ] for
        #                                           inducer_concentration in data]),
        #                            "labels": ["mean",
        #                                       "variance"
        #                                       ]}
        # input_information = {"tf_input": protein_entry,
        #                      "promoter": promoter_entry}
        # promoter_entry = match_promoter(response_characteristic=response_characteristic,
        #                                 reporter_information=reporter_information,
        #                                 input_information=input_information,
        #                                 match_input_sensor=True)

        device_entry = get_device_entry(name=tf_input["group"], group=tf_input["group"], primitive_identifer=["INPUT"],
                                        promoter_id=promoter_entry["identifier"],
                                        # cds_id=protein_entry["identifier"]
                                        )

        if tf_input["group"] == "IPTG":
            reference_sensor_promoter_entry = promoter_entry

        # custom_lib.append(protein_entry)
        custom_lib.append(promoter_entry)
        custom_lib.append(device_entry)
        pass
    return reference_sensor_promoter_entry


def match_promoter(response_characteristic: dict,
                   reporter_information: dict,
                   gate_information: dict = None,
                   input_information: dict = None, match_input_sensor=False):
    """
    This method performs the parameter matching for the promoter model to the response characteristics

    :param response_characteristic: dict providing the input values and corresponding outputcharacteristics
    :return:
    """

    utr_reporter = UTR(reporter_information["cds"]) if "utr" in reporter_information else None
    protein_reporter = Protein(reporter_information["cds"])
    terminator_reporter = None
    reporter_name = protein_reporter.name + "_rpu"
    # Steps
    # Setup Model to train (Gene includes YFP output and the input is a TF Level (possibly derived from RPU))
    # Match model to training data
    # if not match_input_sensor:
    gate_name = None
    if input_information is None:
        raise Exception("Input Information is essentially required for parameter matching.")
    else:
        # tf_input = TranscriptionfactorInputs(input_information["tf_input"])
        promoter_sensor = SensorPromoter(input_information["promoter"])
        input_signal_name = promoter_sensor.cognate_transcription_factors[0]

    genes = []
    if gate_information is None:
        sensor_gene = Gene(id="Sensor", promoter=promoter_sensor, utr=utr_reporter, cds=protein_reporter)
        genes.append(sensor_gene)
        promoter = promoter_sensor
    else:
        gate_name = gate_information["name"]
        utr_gate = UTR(gate_information["utr"]) if "utr" in gate_information else None
        protein_gate = Protein(gate_information["cds"])
        terminator_gate = None
        promoter_gate = Promoter(gate_information["promoter"])

        sensor_gene = Gene(id="Sensor", promoter=promoter_sensor, utr=utr_gate, cds=protein_gate)
        reporter_gene = Gene(id="Reporter", promoter=promoter_gate, utr=utr_reporter, cds=protein_reporter)

        promoter = promoter_gate

        genes.append(sensor_gene)
        genes.append(reporter_gene)

    X = response_characteristic["X"]
    Y = response_characteristic["Y"]
    # X[0] = 10**(-2)
    labels = response_characteristic["labels"]
    y = 0.05
    s = np.max(X) * 0.05
    k_m = k_m = s * y / (1 - y)
    ligand_to_tf_level = lambda x: 1 / (1 + (x / k_m) ** 1)

    # Case distinction
    # Case 1: Parameter estimation of a genetic gate -> both genes are used and matched input sensors and existing reporters are required
    # Case 2: Parameter estimation of an input sensor -> only the reporter gene is required. The input sensor's promoter is directly attached to it.

    def update_model(params):
        new_params = 10 ** params
        new_params[:2] = params[:2]
        promoter.model.insert_params(new_params)

    def loss_func(Y, Y_pred, custom_weights=None):
        weights = np.ones(shape=Y.shape)
        if custom_weights is not None:
            weights = custom_weights
        weights = weights / np.sum(weights)

        loss = 0.0
        # loss += np.sum((weights * np.power(np.log(Y) - np.log(Y_pred), 2))[Y_pred != 0])
        loss += np.sum((weights[:, 0] * np.power(np.log(Y[:, 0]) - np.log(Y_pred[:, 0]), 2)))
        # loss += np.sum((weights[:, 1] * np.power(np.sqrt(Y[:, 1]) - np.sqrt(Y_pred[:, 1]), 2)))

        return loss.item()

    def error_func(params, training_data):

        X, Y = training_data
        update_model(params)

        # X = response_characteristic["X"]
        # Y = response_characteristic["Y"]
        # labels = response_characteristic["labels"]
        Y_pred = [None] * len(X)
        input_vals_dict = {}
        for iX, x in enumerate(X):
            input_vals_dict[input_signal_name] = x
            y_pred = model(plasmid_input_vals=input_vals_dict)
            y_pred = {label: y_pred[iL] for iL, label in enumerate(["mean", "variance"])}
            Y_pred[iX] = [y_pred[label] for label in labels]

        Y_pred = np.array(Y_pred)
        custom_weights = np.ones(shape=Y_pred.shape)  # Equals average weighting.
        custom_weights[0, 0] = 10
        custom_weights[1, 0] = 4
        custom_weights[-2, 0] = 4
        custom_weights[-1, 0] = 10
        custom_weights[:, 1] = 0.1
        custom_weights[0, 1] = 0.2
        custom_weights[-1, 1] = 0.2
        # custom_weights[:, 1] = custom_weights[:, 0] / 10
        # Can realize less weighting of variance results.
        loss = loss_func(Y, Y_pred, custom_weights=custom_weights)
        # deviation = deviation_func(params)
        # loss = deviation
        return loss

    #
    # ref_params = [3.589279542832028, 0.003109276598567825, 479.0328585600103, 5.579281617098493e-11,
    #               5.544387725272706e-11, 5.7219952157354176e-11, 5.581332098763365e-11, 6.669774560347158e-11,
    #               1167.502391260739, 4834.371783185657, 2.1532840435156686e-06, 92474.44318734983]
    #
    # def deviation_func(params):
    #     y_on_ref = ref_params[0]
    #     y_off_ref = ref_params[1]
    #
    #     alphas_ref = ref_params[2:2 + 5]
    #     betas_ref = ref_params[2 + 5:]
    #     y_on, y_off, k_01, k_03, k_10, k_12, k_14, k_21, k_25, k_30, k_34, k_41, k_43, k_45, k_52, k_54 = params
    #     alphas = [
    #         k_03 * k_10 * k_21 * k_41 * k_52 + k_03 * k_10 * k_21 * k_43 * k_52 + k_03 * k_14 * k_21 * k_43 * k_52 + k_03 * k_10 * k_21 * k_41 * k_54 + k_03 * k_10 * k_25 * k_41 * k_54 + k_03 * k_10 * k_21 * k_43 * k_54 + k_03 * k_14 * k_21 * k_43 * k_54 + k_03 * k_10 * k_25 * k_43 * k_54 + k_03 * k_14 * k_25 * k_43 * k_54,
    #         k_01 * k_14 * k_21 * k_30 * k_52 + k_03 * k_10 * k_21 * k_34 * k_52 + k_03 * k_14 * k_21 * k_34 * k_52 + k_01 * k_14 * k_21 * k_43 * k_52 + k_03 * k_10 * k_21 * k_45 * k_52 + k_01 * k_14 * k_21 * k_30 * k_54 + k_01 * k_14 * k_25 * k_30 * k_54 + k_03 * k_10 * k_21 * k_34 * k_54 + k_03 * k_14 * k_21 * k_34 * k_54 + k_03 * k_10 * k_25 * k_34 * k_54 + k_03 * k_14 * k_25 * k_34 * k_54 + k_01 * k_14 * k_21 * k_43 * k_54 + k_03 * k_12 * k_25 * k_43 * k_54 + k_01 * k_14 * k_25 * k_43 * k_54,
    #         k_01 * k_12 * k_25 * k_30 * k_41 + k_03 * k_12 * k_25 * k_34 * k_41 + k_01 * k_12 * k_25 * k_30 * k_43 + k_01 * k_14 * k_21 * k_30 * k_45 + k_01 * k_14 * k_25 * k_30 * k_45 + k_03 * k_10 * k_21 * k_34 * k_45 + k_03 * k_14 * k_21 * k_34 * k_45 + k_03 * k_10 * k_25 * k_34 * k_45 + k_03 * k_14 * k_25 * k_34 * k_45 + k_01 * k_14 * k_21 * k_34 * k_52 + k_01 * k_12 * k_25 * k_30 * k_54 + k_01 * k_14 * k_21 * k_34 * k_54 + k_03 * k_12 * k_25 * k_34 * k_54 + k_01 * k_14 * k_25 * k_34 * k_54 + k_01 * k_12 * k_25 * k_43 * k_54,
    #         k_01 * k_12 * k_25 * k_34 * k_41 + k_01 * k_12 * k_25 * k_30 * k_45 + k_01 * k_14 * k_21 * k_34 * k_45 + k_03 * k_12 * k_25 * k_34 * k_45 + k_01 * k_14 * k_25 * k_34 * k_45 + k_01 * k_12 * k_25 * k_34 * k_54,
    #         k_01 * k_12 * k_25 * k_34 * k_45]
    #     alphas = np.array(alphas)
    #     betas = [
    #         k_10 * k_21 * k_30 * k_41 * k_52 + k_10 * k_21 * k_30 * k_43 * k_52 + k_14 * k_21 * k_30 * k_43 * k_52 + k_10 * k_21 * k_30 * k_41 * k_54 + k_10 * k_25 * k_30 * k_41 * k_54 + k_10 * k_21 * k_30 * k_43 * k_54 + k_14 * k_21 * k_30 * k_43 * k_54 + k_10 * k_25 * k_30 * k_43 * k_54 + k_14 * k_25 * k_30 * k_43 * k_54,
    #         k_01 * k_21 * k_30 * k_41 * k_52 + k_03 * k_21 * k_34 * k_41 * k_52 + k_10 * k_21 * k_34 * k_41 * k_52 + k_01 * k_21 * k_30 * k_43 * k_52 + k_10 * k_21 * k_30 * k_45 * k_52 + k_01 * k_21 * k_30 * k_41 * k_54 + k_01 * k_25 * k_30 * k_41 * k_54 + k_03 * k_21 * k_34 * k_41 * k_54 + k_10 * k_21 * k_34 * k_41 * k_54 + k_03 * k_25 * k_34 * k_41 * k_54 + k_10 * k_25 * k_34 * k_41 * k_54 + k_01 * k_21 * k_30 * k_43 * k_54 + k_01 * k_25 * k_30 * k_43 * k_54 + k_12 * k_25 * k_30 * k_43 * k_54,
    #         k_01 * k_12 * k_30 * k_41 * k_52 + k_03 * k_12 * k_34 * k_41 * k_52 + k_01 * k_21 * k_34 * k_41 * k_52 + k_01 * k_12 * k_30 * k_43 * k_52 + k_01 * k_14 * k_30 * k_45 * k_52 + k_01 * k_21 * k_30 * k_45 * k_52 + k_03 * k_10 * k_34 * k_45 * k_52 + k_03 * k_14 * k_34 * k_45 * k_52 + k_03 * k_21 * k_34 * k_45 * k_52 + k_10 * k_21 * k_34 * k_45 * k_52 + k_01 * k_12 * k_30 * k_41 * k_54 + k_03 * k_12 * k_34 * k_41 * k_54 + k_01 * k_21 * k_34 * k_41 * k_54 + k_01 * k_25 * k_34 * k_41 * k_54 + k_01 * k_12 * k_30 * k_43 * k_54,
    #         k_01 * k_12 * k_34 * k_41 * k_52 + k_01 * k_12 * k_30 * k_45 * k_52 + k_03 * k_12 * k_34 * k_45 * k_52 + k_01 * k_14 * k_34 * k_45 * k_52 + k_01 * k_21 * k_34 * k_45 * k_52 + k_01 * k_12 * k_34 * k_41 * k_54,
    #         k_01 * k_12 * k_34 * k_45 * k_52]
    #     betas = np.array(betas)
    #     divisor = len(alphas) + len(betas)
    #     deviation = np.sum(np.power(alphas - alphas_ref, 2))
    #     deviation += np.sum(np.power(betas - betas_ref, 2))
    #     deviation /= divisor
    #     return deviation

    matching_modes = ["det", "samp"]
    sample_counts = [1, 30]
    matching_modes = ["det"]
    sample_counts = [1]

    best_fun = np.infty
    best_params = None
    for iR in range(5):
        initial_params = None
        for matching_mode, sample_count in zip(matching_modes, sample_counts):
            def model(plasmid_input_vals):
                sim_settings = {"mode": matching_mode}
                output_vals = np.empty(shape=sample_count)
                for iS in range(sample_count):
                    cell_state = {"energy": 0}
                    cell_state.update(plasmid_input_vals)

                    for gene in genes:
                        output_dict = gene(cell_state, sim_settings=sim_settings)
                    output_vals[iS] = cell_state[reporter_name]
                return np.mean(output_vals), np.var(output_vals)  # + 0.1

            # The optimization routine adapts the rates in the infinitesimal generator and the promoter activity

            mask_mean = [label == "mean" for label in labels]
            Y_mean = [y[mask_mean] for y in Y]
            y_on = np.max(Y_mean)
            y_off = np.min(Y_mean)
            # y_on = ref_params[0]
            # y_off = ref_params[1]

            # bounds = [(y_on, y_on * 10), (y_off * 0.1, y_off)] + [[10 ** (-5), 10 ** 5] for _ in range(14)]
            # # Imply inhibitory behaviour via bounds.
            # bounds[3][0] = 2  # k03 > 10
            # bounds[9][1] = 0.5  # k30 < 0.1
            # bounds[8][1] = 0.5  # k25 < 0.1
            # bounds[14][0] = 2  # k52 > 10

            # bounds = [(y_on, y_on), (y_off * 0.5, y_off)] + [[-5, 5] for _ in range(14)]
            bounds = [(y_on, y_on * 2), (y_off * 0.5, y_off)] + [[-5, 5] for _ in range(14)]
            # training_data = [np.array([X[0], X[-1]]), np.array([Y[0], Y[-1]])]
            training_data = [X, Y]

            if initial_params is None:
                initial_params = [y_on, y_off]
                initial_params += list(np.exp(np.random.rand(14) * 2 - 1))
                initial_params = np.array(initial_params)
                new_params = initial_params
                # initial_params[3] = 10
                # initial_params[9] = 0.1
                # initial_params[8] = 0.1
                # initial_params[14] = 10

            if best_params is not None:
                initial_params = best_params
            # loss = error_func(initial_params)
            x0 = initial_params
            result = minimize(error_func, x0=x0,
                              args=training_data,
                              bounds=bounds,
                              method="powell",
                              # method="nelder-mead",
                              # method="TNC",
                              #   method="SLSQP",
                              tol=10 ** (-10),
                              options={"disp": True,
                                       "ftol": 10 ** (-7),
                                       "maxfev": len(x0) * 10 if matching_mode == "samp" else len(x0) * 1500
                                       # len(x0) * 2000
                                       }
                              )
            new_params = result.x
            fun = result.fun
            print(f"Run {iR}: {result.fun} (y_on={y_on}, y_off={y_off})")

            initial_params = new_params

            update_model(new_params)
            # promoter.model.insert_params(new_params)
            matching_mode = "samp"
            sample_count = 100
            Y_pred = [model({input_signal_name: x}) for x in X]
            Y_pred = np.array(Y_pred)

            # plt.figure()
            # ax = plt.gca()
            # ax.plot(X, Y[:, 0], "--k", label="Mean")
            # if Y.shape[1] == 2:
            #     ax.fill_between(X,
            #                     Y[:, 0] - np.sqrt(Y[:, 1]), Y[:, 0] + np.sqrt(Y[:, 1]),
            #                     color="yellow", alpha=0.1)
            #
            # ax.plot(X, Y_pred[:, 0], label="Mean Pred")
            # if Y_pred.shape[1] == 2:
            #     ax.fill_between(X,
            #                     Y_pred[:, 0] - np.sqrt(Y_pred[:, 1]), Y_pred[:, 0] + np.sqrt(Y_pred[:, 1]),
            #                     color="blue", alpha=0.1)
            # ax.set_xscale("log")
            # ax.set_yscale("log")
            # plt.show()

        if fun < best_fun:
            best_fun = fun
            best_params = new_params
    if best_params is not None:
        # promoter.model.insert_params(best_params)
        update_model(new_params)

    print(f"Best Run with error: {best_fun}")
    X_test = np.logspace(-5, np.max(np.log10(X[1:])), 100)
    X_test = X
    # X_tf_level = ligand_to_tf_level(X_test)
    X_tf_level = X
    Y_pred = [model({input_signal_name: x}) for x in X_tf_level]
    Y_pred = np.array(Y_pred)

    plt.figure()
    ax = plt.gca()
    ax.plot(X, Y[:, 0], "--k", label="Mean")
    if Y.shape[1] == 2:
        ax.fill_between(X,
                        Y[:, 0] - np.sqrt(Y[:, 1]), Y[:, 0] + np.sqrt(Y[:, 1]),
                        color="yellow", alpha=0.1)

    ax.plot(X_test, Y_pred[:, 0], label="Mean Pred")
    if Y_pred.shape[1] == 2:
        ax.fill_between(X_test,
                        Y_pred[:, 0] - np.sqrt(Y_pred[:, 1]), Y_pred[:, 0] + np.sqrt(Y_pred[:, 1]),
                        color="blue", alpha=0.1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    dir_path = "ARCTICsim/simulator_nonequilibrium/data/gate_libs/figures/"
    plt.savefig(dir_path + gate_name + ".png", dpi=300)
    # plt.show()

    model_info = promoter_entry["model_info"]
    model_info["PROMOTER_ACTIVITY"] = list(promoter.model.promoter_activity)
    matrices = promoter.model.infinitesimal_generator_function["matrices"]
    model_info["INFINITESIMAL_GENERATOR_FUNCTION"] = {"matrices": {key: matrices[key].tolist() for key in matrices}}
    return promoter_entry


"""
These values are extracted from the paper and the data provided for the figures in the paper.
"""

atc_to_rpu = {0: 0.002,
              5: 0.004,
              10: 0.008,
              15: 0.013,
              20: 0.168,
              22.5: 0.445,
              25: 0.570,
              27.5: 0.668,
              30: 1.019,
              35: 1.227,
              40: 1.616,
              50: 1.963,
              100: 2.455,
              200: 2.548}
atc_to_rpu = {float(key): atc_to_rpu[key] for key in atc_to_rpu}

iptg_to_rpu = {0: 0.011,
               0.1: 0.017,
               0.25: 0.028,
               0.5: 0.049,
               1: 0.132,
               1.5: 0.224,
               2: 0.320,
               2.5: 0.421,
               3: 0.519,
               5: 0.822,
               10: 1.127,
               20: 1.297,
               30: 1.343,
               40: 1.365}
iptg_to_rpu = {float(key): iptg_to_rpu[key] for key in iptg_to_rpu}

xyl_to_rpu = {0: 0.003,
              0.1: 0.006,
              0.25: 0.025,
              0.5: 0.143,
              1: 0.480,
              1.5: 0.859,
              2: 1.133,
              3: 1.409,
              5: 1.662,
              10: 1.804,
              20: 1.904,
              30: 1.949}

xyl_to_rpu = {float(key): xyl_to_rpu[key] for key in xyl_to_rpu}

reference_inducer_inputs = np.array([0, 0.1, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 5, 10, 20])
"""
End of Paper data
"""

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
    cello_toxicity_by_name = {elem["name"]: elem for elem in cello_toxicity}

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
    # Our rate information
    transcription_rate = 24 / 3600
    translation_rate = 200 / 3600
    rna_degradation_rate = 5.8 * 10 ** (-4)
    protein_degradation_rate = 2.9 * 10 ** (-4)
    l_m = translation_rate / rna_degradation_rate
    l_p = translation_rate / protein_degradation_rate
    eu_to_rpu_scaling_factor = (l_m * l_p) ** (-1)


    cello_gates = cello_library_by_collection["gates"]
    custom_lib = []

    # Reporter Output
    rates = {"transcription_rate": transcription_rate,
             "translation_rate": translation_rate,
             "rna_degradation_rate": rna_degradation_rate,
             "protein_degradation_rate": protein_degradation_rate}
    reporters = [{"name": "YFP", "utr_name": "Kozak6"}]
    add_outputs(custom_lib, reporters=reporters, rates=rates)

    # ToDo Extract Report to use for parameter estimation
    reporter_to_use = "YFP"
    reporter_protein = [elem for elem in custom_lib if elem["identifier"] == "protein_" + reporter_to_use][0]
    reporter_information = {"cds": reporter_protein}

    # Transcription Factor Inputs
    tf_inputs = [{"group": "aTc", "promoter": "Ptet", "protein": "tetR", "low": 0.0, "high": 100.0, "data": atc_to_rpu},
                 {"group": "xyl", "promoter": "Pxyl", "protein": "xylR", "low": 0.0, "high": 10.0, "data": xyl_to_rpu},
                 {"group": "IPTG", "promoter": "Plac", "protein": "lacl", "low": 0.0, "high": 20.0, "data": iptg_to_rpu}
                 ]

    reference_sensor_promoter_entry = add_inputs(custom_lib, tf_inputs=tf_inputs,
                                                 reporter_information=reporter_information)

    # ToDo Extract Input to use for parameter estimation

    # Actual gates
    for cello_gate in cello_gates:
        gate_name = cello_gate["name"]
        structure = cello_structures_by_name[cello_gate["structure"]]
        raw_cytometry = cello_cytometry_by_name[gate_name + "_cytometry"]
        toxicity = cello_toxicity_by_name[gate_name + "_toxicity"]
        cytometry_data_dict = GateCytometryData.cello_cytometry_data_to_dict(raw_cytometry[dict_ids["table"]])
        cytometry = GateCytometryData(data_dict=cytometry_data_dict, gate_name=gate_name, cutoff_percentage=0.1)

        promoter_name = structure["outputs"][0].split("_")[0]
        promoter_sequence_names = structure["outputs"]
        utr_name = [elem for elem in structure["devices"] if "cassette" in elem["name"]][0]["components"][0]
        protein_name = cello_gate["regulator"]
        protein_sequence_names = [elem["components"][1] for elem in structure["devices"] if
                                  "cassette" in elem["name"]]

        cello_utr_part = cello_parts_by_name[utr_name]
        cello_promoter_part = [cello_parts_by_name[name] for name in promoter_sequence_names]
        cello_regulator_parts = [cello_parts_by_name[name] for name in protein_sequence_names]

        utr_sequence_entry = get_sequence_entry(name=utr_name, sequence=cello_utr_part["dnasequence"])

        protein_sequence_entries = [get_sequence_entry(name=name, sequence=cello_parts_by_name[name]["dnasequence"])
                                    for name in protein_sequence_names]
        promoter_sequence_entries = [get_sequence_entry(name=name, sequence=cello_parts_by_name[name]["dnasequence"])
                                     for name in promoter_sequence_names]

        utr_entry = {"collection": "utrs",
                     "identifier": "utr_" + utr_name,
                     "name": utr_name,
                     "sequence_id": utr_sequence_entry["identifier"]}

        cds_length = np.median([len(seq["sequence"]) for seq in protein_sequence_entries])
        rna_length = len(utr_sequence_entry["sequence"]) + cds_length
        protein_length = cds_length / 3


        protein_entry = get_protein_entry(name=protein_name,
                                          sequence_ids=[seq["identifier"] for seq in protein_sequence_entries],
                                          transcription_rate=transcription_rate,
                                          rna_degradation_rate=rna_degradation_rate, rna_length=rna_length,
                                          translation_rate=translation_rate,
                                          protein_degradation_rate=protein_degradation_rate,
                                          protein_length=protein_length)

        cognate_transcription_factors = [protein_name, protein_sequence_names]
        promoter_entry = get_promoter_entry(name=promoter_name,
                                            cognate_transcription_factors=cognate_transcription_factors,
                                            sequence_ids=[seq_name for seq_name in promoter_sequence_names],
                                            scaling_factor=eu_to_rpu_scaling_factor)
        # The input levels are in RPU. However, the gene-centric approach requires the input as either inducer levels or TF levels
        # This is achieved by using the inducer concentrations from the paper associated to the respective rpu values.
        # X = cytometry.get_input_levels()
        X = reference_inducer_inputs
        Y = [cytometry.get_histogram(x) for x in cytometry.get_input_levels()]
        Y = [(hist.mean(), hist.variance()) for hist in Y]
        response_characteristic = {"X": np.array(X),
                                   "Y": np.array(Y),
                                   "labels": ["mean", "variance"]}
        input_information = {"promoter": reference_sensor_promoter_entry}
        gate_information = {  # "utr":None,
            "cds": protein_entry,
            "promoter": promoter_entry,
            "name": gate_name}
        promoter_entry = match_promoter(response_characteristic=response_characteristic,
                                        reporter_information=reporter_information,
                                        input_information=input_information,
                                        gate_information=gate_information,
                                        match_input_sensor=False)

        device_entry = get_device_entry(name=gate_name, group=cello_gate["group"], regulator=cello_gate["regulator"],
                                        primitive_identifer=["NOT", "NOR2"], utr_id=utr_entry["identifier"],
                                        cds_id=protein_entry["identifier"], promoter_id=promoter_entry["identifier"],
                                        color=cello_gate["color"])

        custom_lib += [device_entry, utr_entry, protein_entry, promoter_entry]
        # custom_lib += [device_entry, protein_entry, promoter_entry]

        # Next steps:
        # Match Characteristics of promoter to cytometry data
        # Infer Inflection Point
        # Store everything in a list

        # Add everything to the lists

        pass
        pass
        # Write out everything to a json file


    custom_lib = {elem["identifier"]: elem for elem in custom_lib}
    custom_lib = list(custom_lib.values())

    output_dir = "ARCTICsim/simulator_nonequilibrium/data/gate_libs/"
    with open(output_dir + "gate_lib_yeast_generated.json", "w") as file:
        json.dump(custom_lib, file, indent=4)
