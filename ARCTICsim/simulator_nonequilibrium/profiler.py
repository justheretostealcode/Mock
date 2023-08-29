import json
import cProfile

from datetime import datetime

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
    file_name= f"profile_results_{formatted_timestamp}.prof"
    profiler.dump_stats(file_name)
    print(f"Visualize results with: snakeviz {file_name}")
    pass


def profile_four_state_promoter_model(n=100):
    gate_lib_path = "data/gate_libs/gate_lib_yeast.json"

    with open(gate_lib_path, "r") as file:
        gate_lib_data = json.load(file)

    not_gates = [gate for gate in gate_lib_data if "NOT" in gate["primitiveIdentifier"]]

    test_gate_entry = not_gates[0]

    model_entry = test_gate_entry["biorep"]["model"]

    model = FourStatePromoterModel(model_entry)

    profiler = start_profiler()

    for iN in range(n):
        model.get_distributions_and_energy_rate(10000)

    stop_profiler(profiler)


if __name__ == '__main__':
    n = 1000

    profile_four_state_promoter_model(n=n)
