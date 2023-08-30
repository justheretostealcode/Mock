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

    results = None
    results_file_path = "results.json"
    if os.path.exists(results_file_path):
        with open(results_file_path, "r") as file:
            results = json.load(file)
        results = np.array(results)

    not_gates = [gate for gate in gate_lib_data if "NOT" in gate["primitiveIdentifier"]]

    test_gate_entry = not_gates[0]

    model_entry = test_gate_entry["biorep"]["model"]

    model = FourStatePromoterModel(model_entry)

    profiler = start_profiler()

    compare_results = results is not None

    cur_results = np.empty(shape=(n, 5))

    for iN in range(n):
        result = model.get_distributions_and_energy_rate(10 ** 2 * iN * np.pi + 123.123451312431 * np.pi)

        cur_results[iN] = result

    stop_profiler(profiler)

    if compare_results:
        diffs = cur_results - results
        print("Diff Error:", np.mean(np.square(diffs)))
    else:
        results = list(map(list, cur_results))
        with open(results_file_path, "w") as file:
            json.dump(results, file)


if __name__ == '__main__':
    n = 10000

    profile_four_state_promoter_model(n=n)
