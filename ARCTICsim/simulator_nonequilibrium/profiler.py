import json
import cProfile
import os

from datetime import datetime

import numpy as np

from models.four_state_promoter_model import FourStatePromoterModel


def start_profiler():
    profiler = cProfile.Profile()
    profiler.enable()
    return profiler


def stop_profiler(profiler):
    # Get the current timestamp
    current_timestamp = datetime.now()
    formatted_timestamp = current_timestamp.strftime("%Y-%m-%d_%H-%M-%S")
    profiler.disable()
    file_name = f"profile_results_{formatted_timestamp}.prof"
    profiler.dump_stats(file_name)
    print(f"Visualize results with: snakeviz {file_name}")
    pass


def profile_four_state_promoter_model(n=100):
    gate_lib_path = "data/gate_libs/gate_lib_yeast.json"

    with open(gate_lib_path, "r") as file:
        gate_lib_data = json.load(file)

    results1 = None
    results2 = None
    results_file_path = "results.json"
    if os.path.exists(results_file_path):
        with open(results_file_path, "r") as file:
            input_data = json.load(file)


        results1 = np.array(input_data[0])
        results2 = np.array(input_data[1])

    not_gates = [gate for gate in gate_lib_data if "NOT" in gate["primitiveIdentifier"]]



    model_entry_1 = not_gates[0]["biorep"]["model"]
    model_entry_2 = not_gates[1]["biorep"]["model"]

    model_1 = FourStatePromoterModel(model_entry_1)
    model_2 = FourStatePromoterModel(model_entry_2)

    profiler = start_profiler()

    compare_results = results1 is not None and results2 is not None

    cur_results1 = np.empty(shape=(n, 5))
    cur_results2 = np.empty(shape=(n, 5))

    for iN in range(n):
        in_val = 10 ** 3 * iN * np.pi + 123.123451312431 * np.pi
        result1 = model_1.get_distributions_and_energy_rate(in_val)
        result2 = model_2.get_distributions_and_energy_rate(in_val)

        cur_results1[iN] = result1
        cur_results2[iN] = result2

    stop_profiler(profiler)

    if compare_results:
        diffs1 = cur_results1 - results1
        diffs2 = cur_results2 - results2
        print("Diff Error 1:", np.mean(np.square(diffs1)))
        print("Diff Error 2:", np.mean(np.square(diffs2)))
    else:
        cur_results1 = list(map(list, cur_results1))
        cur_results2 = list(map(list, cur_results2))
        output_data = [cur_results1, cur_results2]
        with open(results_file_path, "w") as file:
            json.dump(output_data, file)


if __name__ == '__main__':
    n = 10000

    profile_four_state_promoter_model(n=n)
