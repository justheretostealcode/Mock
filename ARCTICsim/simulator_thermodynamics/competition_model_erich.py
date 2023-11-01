import json
import matplotlib.pyplot as plt

# Threshold to be used to check for maximal allowed resource demand
THRESHOLD = 0.8  # TODO


def add_gate_resource_competition(
    library: str, Q_max: float, Q_min: float, S: float, name: str
) -> None:
    """
    This function adds resource competition data of a gate in lib file.

    Parameters:
    library (json file path): Json file path of library.
    Q_max (float): maximum demand of gate.
    Q_min (float): minimum demand of gate.
    S (float): Sensitivity of gate.
    name (str): name of gate to modify

    Returns:
    None: The lib was successfully modified.
    """
    with open(library, "r") as file:
        lib = json.load(file)
    lib_gates = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "gates"
    ][0]
    gate_to_modify = [gate for gate in lib_gates if gate["name"] == name][0]
    gate_to_modify["parameters"]["Q_max"] = Q_max
    gate_to_modify["parameters"]["Q_min"] = Q_min
    gate_to_modify["parameters"]["S"] = S
    with open(library, "w") as file:
        json.dump(lib, file)


def add_output_resource_competition(
    library: str, Q_y: float, S: float, name: str
) -> None:
    """
    This function adds resource competition data of an output in lib file.

    Parameters:
    library (json file path): Json file path of library.
    Q_y (float):  demand slope of gate.
    S (float): Sensitivity of gate.
    name (str): name of gate to modify

    Returns:
    None: The lib was successfully modified.
    """
    with open(library, "r") as file:
        lib = json.load(file)
    lib_gates = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "gates"
    ][0]
    gate_to_modify = [gate for gate in lib_gates if gate["name"] == name][0]
    if not gate_to_modify["associated_devices"][0].startswith("output"):
        raise ValueError("Not an output gate!")
    gate_to_modify["parameters"]["Q_y"] = Q_y
    gate_to_modify["parameters"]["S"] = S
    with open(library, "w") as file:
        json.dump(lib, file)


def add_input_resource_competition(library: str, Q: float, S: float, name: str) -> None:
    """
    This function adds resource competition data of an output in lib file.

    Parameters:
    library (json file path): Json file path of library.
    Q (float):  demand of gate.
    S (float): Sensitivity of gate.
    name (str): name of gate to modify

    Returns:
    None: The lib was successfully modified.
    """
    with open(library, "r") as file:
        lib = json.load(file)
    lib_inputs = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "transcription_factors"
    ][0]
    lib_inputs = [
        inp for inp in lib_inputs if inp["associated_devices"][0].startswith("input")
    ]
    gate_to_modify = [gate for gate in lib_inputs if gate["name"] == name][0]
    gate_to_modify["parameters"] = {}
    gate_to_modify["parameters"]["Q"] = Q
    gate_to_modify["parameters"]["S"] = S
    with open(library, "w") as file:
        json.dump(lib, file)


def _add_test_data(library: str) -> None:
    """
    This function adds hill function data in lib file.

    Parameters:
    library (json file path): Json file path of library.

    Returns:
    None: The lib was successfully modified.
    """
    with open(library, "r") as file:
        lib = json.load(file)
    lib_gates = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "gates"
    ][0]
    lib_inputs = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "transcription_factors"
    ][0]
    lib_inputs = [
        inp for inp in lib_inputs if inp["associated_devices"][0].startswith("input")
    ]
    # delete old parameters and add test parameters
    for gate in lib_gates:
        if gate["associated_devices"][0].startswith("output"):
            del gate["parameters"]["mu"]
            del gate["parameters"]["gamma"]
            del gate["parameters"]["delta"]
            gate["parameters"]["Q_y"] = 0.8
            gate["parameters"]["S"] = 0.5
        else:
            del gate["parameters"]["mu"]
            del gate["parameters"]["gamma"]
            del gate["parameters"]["delta"]
            del gate["parameters"]["beta"]
            del gate["parameters"]["mu_p"]
            del gate["parameters"]["gamma_p"]
            del gate["parameters"]["kappa"]
            del gate["parameters"]["n"]
            gate["parameters"]["y_max"] = 4
            gate["parameters"]["y_min"] = 0.01
            gate["parameters"]["K"] = 0.1
            gate["parameters"]["n"] = 3
            gate["parameters"]["Q_max"] = 0.2
            gate["parameters"]["Q_min"] = 0.05
            gate["parameters"]["S"] = 0.1
    for inp in lib_inputs:
        inp["parameters"] = {"Q": 0.05, "S": 0.01}
    # with open(
    #     "./../ARCTICsim/erich_libs/gate_lib_demand_sensitivity.json", "w"
    # ) as json_file:
    #     json.dump(lib, json_file)


# _add_test_data("./../ARCTICsim/iwbda_libs/iwbda_lib_id_binary_delta.json")


def _add_test_data_for_simulation(library: str) -> None:
    """
    This function adds hill function data in lib file.

    Parameters:
    library (json file path): Json file path of library.

    Returns:
    None: The lib was successfully modified.
    """
    with open(library, "r") as file:
        lib = json.load(file)
    lib_gates = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "gates"
    ][0]
    lib_inputs = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "transcription_factors"
    ][0]
    lib_inputs = [
        inp for inp in lib_inputs if inp["associated_devices"][0].startswith("input")
    ]
    # delete old parameters and add test parameters
    for gate in lib_gates:
        if gate["associated_devices"][0].startswith("output"):
            del gate["parameters"]["mu"]
            del gate["parameters"]["gamma"]
            del gate["parameters"]["delta"]
            gate["parameters"]["Q_y"] = 0.8
            gate["parameters"]["S"] = 0.5
            gate["levels"] = {"on": 0}
            gate["levels"]["off"] = 0
        else:
            del gate["parameters"]["mu"]
            del gate["parameters"]["gamma"]
            del gate["parameters"]["delta"]
            del gate["parameters"]["beta"]
            del gate["parameters"]["mu_p"]
            del gate["parameters"]["gamma_p"]
            del gate["parameters"]["kappa"]
            del gate["parameters"]["n"]
            gate["levels"] = {"on": 4}
            gate["levels"]["off"] = 0.01
            gate["hill_parameters"] = {"K": 0.1}
            gate["hill_parameters"]["n"] = 3
            gate["parameters"]["Q_max"] = 0.2
            gate["parameters"]["Q_min"] = 0.05
            gate["parameters"]["S"] = 0.1
    for inp in lib_inputs:

        inp["parameters"] = {"Q": 0.05, "S": 0.01}
    with open(
        "./../ARCTICsim/erich_libs/gate_lib_demand_sensitivity.json", "w"
    ) as json_file:
        json.dump(lib, json_file)


# _add_test_data_for_simulation("./../ARCTICsim/iwbda_libs/iwbda_lib_id_binary_delta.json")


def _preprocess_circuit_json(assigns: dict, lib_gates: list, struct: dict) -> dict:
    """
    This function combines in a preprocess step the json files of the simulation.

    Parameters:
    assigns (dict): dict of assignments.
    lib_gates (list): list of library gates.
    struct (dict): dict of circuit topology.

    Returns:
    dict: combined dict of the files.
    """
    edges = [edge for edge in struct["graph"]["edges"]]
    predecessors = {
        logic_target["target"]: {
            "sources": [
                edge["source"]
                for edge in edges
                if edge["target"] == logic_target["target"]
            ],
            "output": [
                truth_values
                for logic_gate, truth_values in struct["gate_truthtables"].items()
                if logic_gate == logic_target["target"]
            ][0],
        }
        for logic_target in edges
    }
    # predecessors.pop("O")
    for target, values in predecessors.items():
        for logic_gate, assigned_gate in assigns.items():
            if logic_gate == target:
                predecessors[target]["assigned"] = assigns[logic_gate]
    for target, values in predecessors.items():
        for lib_gate in lib_gates:
            if lib_gate["associated_devices"][0] == values["assigned"]:
                if lib_gate["associated_devices"][0].startswith("output"):
                    predecessors[target]["parameters"] = {
                        "Q_y": lib_gate["parameters"]["Q_y"],
                        "S": lib_gate["parameters"]["S"],
                    }
                else:
                    predecessors[target]["parameters"] = {
                        "Q_max": lib_gate["parameters"]["Q_max"],
                        "Q_min": lib_gate["parameters"]["Q_min"],
                        "S": lib_gate["parameters"]["S"],
                        "n": lib_gate["parameters"]["n"],
                        "K": lib_gate["parameters"]["K"],
                        "y_max": lib_gate["parameters"]["y_max"],
                        "y_min": lib_gate["parameters"]["y_min"],
                    }
    [values.pop("assigned") for target, values in predecessors.items()]
    return predecessors


# with open(
#     "./../ARCTICsim/test_cases/bounding-debug/cello_00001110_assignment_optimal.json",
#     "r",
# ) as file:
#     assigns = json.load(file)
# with open("./../ARCTICsim/erich_libs/gate_lib_demand_sensitivity.json", "r") as file:
#     lib = json.load(file)
# with open("./../ARCTICsim/test_cases/bounding-debug/cello_00001110.json", "r") as file:
#     struct = json.load(file)
# lib_gates = [
#     [gate for gate in section["members"]]
#     for section in lib
#     if section["class"] == "gates"
# ][0]
# predecessors = _preprocess_circuit_json(assigns, lib_gates, struct)
# print(predecessors)


def _circuit_demand_estimate(assignments: str, library: str, structure: str) -> bool:
    """
    This function estimates the circuit demand and compares it to the threshold.

    Parameters:
    assignments (json file path): json file path of assignments.
    library (json file path): json file path of library.
    structure (json file path): json file path of structure.

    Returns:
    bool: True, if circuit is accepted, else False.
    """
    with open(assignments, "r") as file:
        assigns = json.load(file)
    with open(library, "r") as file:
        lib = json.load(file)
    with open(structure, "r") as file:
        struct = json.load(file)
    lib_gates = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "gates"
    ][0]
    predecessors = _preprocess_circuit_json(assigns, lib_gates, struct)
    # sum up demand of all gates (without output/input)
    Q_max = [
        sum(
            [
                values["parameters"]["Q_max"]
                if values["output"][i] == "0"
                else (
                    2 * values["parameters"]["Q_min"]
                    if len(values["sources"]) == 2
                    else values["parameters"]["Q_min"]
                )
            ][0]
            for target, values in predecessors.items()
            if target != "O"
        )
        for i in range(len(next(iter(struct["gate_truthtables"].values()))))
    ]
    # add output resource demand
    for i in range(len(next(iter(struct["gate_truthtables"].values())))):
        for target, values in predecessors.items():
            if target == "O":
                if values["sources"][0] in ["a", "b", "c"]:  # TODO
                    input1 = 1
                else:
                    for gateName, gateVals in predecessors.items():
                        if gateName == values["sources"][0]:
                            if gateVals["output"][i] == "1":
                                input1 = gateVals["parameters"]["y_max"]
                            else:
                                input1 = gateVals["parameters"]["y_min"]
                if len(values["sources"]) == 2:
                    if values["sources"][1] in ["a", "b", "c"]:  # TODO
                        input2 = 0
                    else:
                        for gateName, gateVals in predecessors.items():
                            if gateName == values["sources"][1]:
                                if gateVals["output"][i] == "1":
                                    input2 = gateVals["parameters"]["y_max"]
                                else:
                                    input2 = gateVals["parameters"]["y_min"]
                else:
                    input2 = 0
                Q_max[i] += values["parameters"]["Q_y"] * (input1 + input2)
    # add constant demand of constitutive input promoters of the circuit. (tac,bad,tet)
    # Annahme: jeder Circuit hat immer diese inputs a,b,c.
    lib_inputs = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "transcription_factors"
    ][0]
    lib_inputs = [
        inp for inp in lib_inputs if inp["associated_devices"][0].startswith("input")
    ]
    for i in range(len(next(iter(struct["gate_truthtables"].values())))):
        Q_max[i] += sum([gate["parameters"]["Q"] for gate in lib_inputs])
    return max(Q_max) < THRESHOLD


# _circuit_demand_estimate(
#     "./../ARCTICsim/test_cases/bounding-debug/cello_00001110_assignment_optimal.json",
#     "./../ARCTICsim/erich_libs/gate_lib_demand_sensitivity.json",
#     "./../ARCTICsim/test_cases/bounding-debug/cello_00001110.json"
# )


def _get_circuit_demand_estimate(
    assigns, lib_gates, struct, lib_inputs, i: int
) -> float:
    """
    This function estimates the circuit demand and compares it to the threshold.

    Parameters:
    assignments (read json file): json file of assignments.
    lib_gates (read json file): json file of gate and output library.
    lib_inputs (read json file): json file of inputs library.
    structure (read json file): json file of structure.
    i (int): i-th input combination.

    Returns:
    float: circuit demand.
    """
    predecessors = _preprocess_circuit_json(assigns, lib_gates, struct)
    # sum up demand of all gates (without output/input)
    Q_c = [
        [
            values["parameters"]["Q_max"]
            if values["output"][i] == "0"
            else (
                2 * values["parameters"]["Q_min"]
                if len(values["sources"]) == 2
                else values["parameters"]["Q_min"]
            )
        ][0]
        for target, values in predecessors.items()
        if target != "O"
    ][0]
    # add output resource demand
    for target, values in predecessors.items():
        if target == "O":
            if values["sources"][0] in ["a", "b", "c"]:  # TODO
                input1 = 1
            else:
                for gateName, gateVals in predecessors.items():
                    if gateName == values["sources"][0]:
                        if gateVals["output"][i] == "1":
                            input1 = gateVals["parameters"]["y_max"]
                        else:
                            input1 = gateVals["parameters"]["y_min"]
            if len(values["sources"]) == 2:
                if values["sources"][1] in ["a", "b", "c"]:  # TODO
                    input2 = 0
                else:
                    for gateName, gateVals in predecessors.items():
                        if gateName == values["sources"][1]:
                            if gateVals["output"][i] == "1":
                                input2 = gateVals["parameters"]["y_max"]
                            else:
                                input2 = gateVals["parameters"]["y_min"]
            else:
                input2 = 0
            Q_c += values["parameters"]["Q_y"] * (input1 + input2)
    # add constant demand of constitutive input promoters of the circuit. (tac,bad,tet)
    # Annahme: jeder Circuit hat immer diese inputs a,b,c.
    Q_c += sum([gate["parameters"]["Q"] for gate in lib_inputs])
    return Q_c


def _inverse_hill_fct(
    y: float, y_max: float, y_min: float, K: float, n: float
) -> float:
    """
    This function computes inverse hill function.

    Parameters:
    y (float): function value to which the x value of the hill function should be computed.
    y_max (float): maximum output of gate.
    y_min (float): minimum output of gate.
    K (float): Threshold of gate.
    n (float): hill coefficient of gate.

    Returns:
    float: x value of hill fct.
    """
    return K * (((y_max - y_min) / (y - y_min) - 1) ** (1.0 / n))


def _max_circuit_sensitivity(
    assigns: dict,
    lib_gates: list,
    struct: dict,
    i: int,
    input_a: list,
    input_b: list,
    input_c: list,
    max_output_0: float,
    min_output_1: float,
    lib_inputs: list,
) -> bool:
    """
    This function fulfills an input threshold analysis with perturbated inputs.

    Parameters:
    assigns (dict): dict of assignments.
    lib_gates (list): list of library gates.
    struct (dict): dict of structure.
    i (int): input combination iteration state index.
    input_a (list): [y_min, y_max] of input 'a'.
    input_b (list): [y_min, y_max] of input 'b'.
    input_c (list): [y_min, y_max] of input 'c'.
    max_output_0 (float): maximum allowed output to be recognized as 0.
    min_output_1 (float): minimum needed output to be recognized as 1.
    lib_inputs (list): list of circuit lib inputs.

    Returns:
    bool: True, if circuit is accepted, else False.
    """
    # Notation: input_x=[y_min_x, y_max_x]
    predecessors = _preprocess_circuit_json(assigns, lib_gates, struct)
    # Note: Q_c muss noch um das jeweilige Q_M reduziert werden. Für jeden output.
    Q_c = _get_circuit_demand_estimate(assigns, lib_gates, struct, lib_inputs, i)
    output = predecessors.pop("O")
    # perturbate inputs
    for lib_inp in lib_inputs:
        for assgn_inp, associated_inp in assigns.items():
            if associated_inp == lib_inp["associated_devices"][0]:
                if assgn_inp == "a":
                    Q_P = Q_c - lib_inp["parameters"]["Q"]
                    input_a[0] = (
                        input_a[0]
                        * (1 + lib_inp["parameters"]["S"] * Q_P)
                        * (1 + lib_inp["parameters"]["Q"])
                        / (1 + lib_inp["parameters"]["Q"] + Q_P)
                    )
                    input_a[1] = (
                        input_a[1]
                        * (1 + lib_inp["parameters"]["S"] * Q_P)
                        * (1 + lib_inp["parameters"]["Q"])
                        / (1 + lib_inp["parameters"]["Q"] + Q_P)
                    )
                elif assgn_inp == "b":
                    Q_P = Q_c - lib_inp["parameters"]["Q"]
                    input_b[0] = (
                        input_b[0]
                        * (1 + lib_inp["parameters"]["S"] * Q_P)
                        * (1 + lib_inp["parameters"]["Q"])
                        / (1 + lib_inp["parameters"]["Q"] + Q_P)
                    )
                    input_b[1] = (
                        input_b[1]
                        * (1 + lib_inp["parameters"]["S"] * Q_P)
                        * (1 + lib_inp["parameters"]["Q"])
                        / (1 + lib_inp["parameters"]["Q"] + Q_P)
                    )
                else:
                    Q_P = Q_c - lib_inp["parameters"]["Q"]
                    input_c[0] = (
                        input_c[0]
                        * (1 + lib_inp["parameters"]["S"] * Q_P)
                        * (1 + lib_inp["parameters"]["Q"])
                        / (1 + lib_inp["parameters"]["Q"] + Q_P)
                    )
                    input_c[1] = (
                        input_c[1]
                        * (1 + lib_inp["parameters"]["S"] * Q_P)
                        * (1 + lib_inp["parameters"]["Q"])
                        / (1 + lib_inp["parameters"]["Q"] + Q_P)
                    )
    # compute and perturbate circuit
    for target, values in predecessors.items():
        inp_1 = values["sources"][0]
        state = int(
            [
                truth_values[i]
                for logic_gate, truth_values in struct["gate_truthtables"].items()
                if logic_gate == inp_1
            ][0]
        )
        if inp_1 == "a":
            inp_1 = input_a[state]
        elif inp_1 == "b":
            inp_1 = input_b[state]
        elif inp_1 == "c":
            inp_1 = input_c[state]
        else:
            Q_M = [
                [
                    values["parameters"]["Q_max"]
                    if state == 0
                    else (
                        2 * values["parameters"]["Q_min"]
                        if len(values["sources"]) == 2
                        else values["parameters"]["Q_min"]
                    )
                ][0]
                for target, values in predecessors.items()
                if target == inp_1
            ][0]
            Q_P = Q_c - Q_M
            S_M = [
                values["parameters"]["S"]
                for target, values in predecessors.items()
                if target == inp_1
            ][0]
            y = [
                [
                    values["parameters"]["y_min"]
                    if state == 0
                    else values["parameters"]["y_max"]
                ][0]
                for target, values in predecessors.items()
                if target == inp_1
            ][0]
            inp_1 = y * (1 + S_M * Q_P) * (1 + Q_M) / (1 + Q_M + Q_P)
        if len(values["sources"]) == 2:
            inp_2 = values["sources"][1]
            state = int(
                [
                    truth_values[i]
                    for logic_gate, truth_values in struct["gate_truthtables"].items()
                    if logic_gate == inp_2
                ][0]
            )
            if inp_2 == "a":
                inp_2 = input_a[state]
            elif inp_2 == "b":
                inp_2 = input_b[state]
            elif inp_2 == "c":
                inp_2 = input_c[state]
            else:
                Q_M = [
                    [
                        values["parameters"]["Q_max"]
                        if state == 0
                        else values["parameters"]["Q_min"]
                    ][0]
                    for target, values in predecessors.items()
                    if target == inp_2
                ][0]
                Q_P = Q_c - Q_M
                S_M = [
                    values["parameters"]["S"]
                    for target, values in predecessors.items()
                    if target == inp_2
                ][0]
                y = [
                    [
                        values["parameters"]["y_min"]
                        if state == 0
                        else values["parameters"]["y_max"]
                    ][0]
                    for target, values in predecessors.items()
                    if target == inp_2
                ][0]
                inp_2 = y * (1 + S_M * Q_P) * (1 + Q_M) / (1 + Q_M + Q_P)
        else:
            inp_2 = 0
        inp = inp_1 + inp_2
        if values["output"][i] == "1":
            max_allowed_y = _inverse_hill_fct(
                0.5 * values["parameters"]["y_max"],
                values["parameters"]["y_max"],
                values["parameters"]["y_min"],
                values["parameters"]["K"],
                values["parameters"]["n"],
            )
            if inp > max_allowed_y:
                return False
        else:
            min_needed_y = _inverse_hill_fct(
                2 * values["parameters"]["y_min"],
                values["parameters"]["y_max"],
                values["parameters"]["y_min"],
                values["parameters"]["K"],
                values["parameters"]["n"],
            )
            if inp < min_needed_y:
                return False
    # check for desired output
    if output["sources"][0] in ["a", "b", "c"]:  # TODO
        input1 = 1
    else:
        for gateName, gateVals in predecessors.items():
            if gateName == output["sources"][0]:
                if gateVals["output"][i] == "1":
                    input1 = gateVals["parameters"]["y_max"]
                else:
                    input1 = gateVals["parameters"]["y_min"]
    if len(output["sources"]) == 2:
        if output["sources"][1] in ["a", "b", "c"]:  # TODO
            input2 = 0
        else:
            for gateName, gateVals in predecessors.items():
                if gateName == output["sources"][1]:
                    if gateVals["output"][i] == "1":
                        input2 = gateVals["parameters"]["y_max"]
                    else:
                        input2 = gateVals["parameters"]["y_min"]
    else:
        input2 = 0
    input_sum = input1 + input2
    Q_P = Q_c - output["parameters"]["Q_y"] * input_sum
    y_M_P = (
        input_sum
        * (1 + output["parameters"]["S"] * Q_P)
        * (1 + output["parameters"]["Q_y"] * input_sum)
        / (1 + output["parameters"]["Q_y"] * input_sum + Q_P)
    )
    if output["output"][i] == "0":
        if y_M_P > max_output_0:
            return False
    else:
        if y_M_P < min_output_1:
            return False
    return True


def _max_circuit_sensitivity_wrapper(
    assignments: str,
    library: str,
    structure: str,
    input_a: list,
    input_b: list,
    input_c: list,
    max_output_0: float,
    min_output_1: float,
) -> bool:
    """
    This function fulfills an input threshold analysis with perturbated inputs.

    Parameters:
    assignments (json file path): json file path of assignments.
    library (json file path): json file path of library.
    structure (json file path): json file path of structure.
    input_a (list): [y_min, y_max] of input 'a'.
    input_b (list): [y_min, y_max] of input 'b'.
    input_c (list): [y_min, y_max] of input 'c'.
    max_output_0 (float): maximum allowed output to be recognized as 0.
    min_output_1 (float): minimum needed output to be recognized as 1.

    Returns:
    bool: True, if circuit is accepted, else False.
    """
    with open(assignments, "r") as file:
        assigns = json.load(file)
    with open(library, "r") as file:
        lib = json.load(file)
    with open(structure, "r") as file:
        struct = json.load(file)
    lib_gates = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "gates"
    ][0]
    # Annahme: jeder Circuit hat immer diese inputs a,b,c.
    lib_inputs = [
        [gate for gate in section["members"]]
        for section in lib
        if section["class"] == "transcription_factors"
    ][0]
    lib_inputs = [
        inp for inp in lib_inputs if inp["associated_devices"][0].startswith("input")
    ]
    for i in range(len(next(iter(struct["gate_truthtables"].values())))):
        if not _max_circuit_sensitivity(
            assigns,
            lib_gates,
            struct,
            i,
            input_a.copy(),
            input_b.copy(),
            input_c.copy(),
            max_output_0,
            min_output_1,
            lib_inputs,
        ):
            return False
    return True


# _max_circuit_sensitivity_wrapper(
#     "./../ARCTICsim/test_cases/bounding-debug/cello_00001110_assignment_optimal.json",
#     "./../ARCTICsim/erich_libs/gate_lib_demand_sensitivity.json",
#     "./../ARCTICsim/test_cases/bounding-debug/cello_00001110.json",
#     [0.01, 4],
#     [0.01, 4],
#     [0.01, 4],
#     0.1,
#     3,
# )


def check_resource_competition(
    assignments: str,
    library: str,
    structure: str,
    input_a: list,
    input_b: list,
    input_c: list,
    max_output_0: float,
    min_output_1: float,
) -> bool:
    """
    This function checks if the circuits resource competition characteristics are okay.

    Parameters:
    assignments (json file path): json file path of assignments.
    library (json file path): json file path of library.
    structure (json file path): json file path of structure.
    input_a (list): [y_min, y_max] of input 'a'.
    input_b (list): [y_min, y_max] of input 'b'.
    input_c (list): [y_min, y_max] of input 'c'.
    max_output_0 (float): maximum allowed output to be recognized as 0.
    min_output_1 (float): minimum needed output to be recognized as 1.

    Returns:
    bool: True, if circuit is accepted, else False.
    """
    return _max_circuit_sensitivity_wrapper(
        assignments,
        library,
        structure,
        input_a,
        input_b,
        input_c,
        max_output_0,
        min_output_1,
    ) and _circuit_demand_estimate(assignments, library, structure)


# check_resource_competition(
#     "./../ARCTICsim/test_cases/bounding-debug/cello_00001110_assignment_optimal.json",
#     "./../ARCTICsim/erich_libs/gate_lib_demand_sensitivity.json",
#     "./../ARCTICsim/test_cases/bounding-debug/cello_00001110.json",
#     [0.01, 4],
#     [0.01, 4],
#     [0.01, 4],
#     0.1,
#     3,
# )


def get_Q_output(y: float, Q_y: float) -> float:
    return Q_y * y


def get_Q_gate(
    y: float, y_max: float, y_min: float, Q_max: float, Q_min: float
) -> float:
    return Q_max - (Q_max - Q_min) * (y - y_min) / (y_max - y_min)


def get_S(y: float, y_0: float, Q_M, Q_p: float) -> float:
    return 0.1
    # return 0 if Q_p == 0 else (y / y_0 * (1 + Q_p + Q_M) / (1 + Q_M) - 1) / Q_p


def perturbed_gate_output_of_gate(
    y: float,
    Q_p: float,
    S: float,
    y_max: float,
    y_min: float,
    Q_max: float,
    Q_min: float,
    y_0: float,
) -> list:
    """
    This function computes the perturbed gate output out of the simulated output
    using resource parameters in the sense of a fixed-point iteration.

    Parameters:
    y (float): simulated output of gate of that iteration.
    Q_p (float): perturbing gate demand of that iteration.
    S (float): gate sensitivity of that iteration.
    y_max (float): max gate output.
    y_min (float): min gate output.
    n (float): hill coefficient of gate.
    K (float): hill param.
    Q_max (float): max demand of gate.
    Q_min (float): min demand of gate.
    y_0 (dict): unperturbed simulation output of gate.

    Returns:
    list: [y_p, Q_M, S]
    """
    Q_M = get_Q_gate(y, y_max, y_min, Q_max, Q_min)
    y_p = y * (1 + S * Q_p) * (1 + Q_M) / (1 + Q_M + Q_p)
    Q_M_next = get_Q_gate(y_p, y_max, y_min, Q_max, Q_min)
    return [y_p, Q_M_next, get_S(y_p, y_0, Q_M_next, Q_p)]


def perturbed_gate_output_of_output(
    y: float, Q_p: float, S: float, Q_y: float, y_0: float
) -> list:
    """
    This function computes the perturbed gate output out of the simulated output
    using resource parameters in the sense of a fixed-point iteration.

    Parameters:
    y (float): simulated output of gate of that iteration.
    Q_p (float): perturbing gate demand of that iteration.
    S (float): gate sensitivity of that iteration.
    Q_y (float): demand per output of gate. #slope!
    y_0 (dict): unperturbed simulation output of gate.

    Returns:
    list: [y_p, Q_M, S]
    """
    Q_M = get_Q_output(y, Q_y)
    y_p = y * (1 + S * Q_p) * (1 + Q_M) / (1 + Q_M + Q_p)
    Q_M_next = get_Q_output(y_p, Q_y)
    return [y_p, Q_M_next, get_S(y_p, y_0, Q_M_next, Q_p)]
