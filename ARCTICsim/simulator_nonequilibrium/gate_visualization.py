import matplotlib
import numpy as np

from simulator.gatelib import GateLib
from simulator.utils import JsonFile
from models.four_state_promoter_model import FourStatePromoterModel

import scipy.stats as sts
import matplotlib.pyplot as plt
import seaborn as sns


font = {
    'size': 14
}
matplotlib.rc('font', **font)

def visualize_protein_distribution(model):
    def log_normal_params(empiric_mean, empiric_var):
        common_val = empiric_var / empiric_mean ** 2 + 1
        mu = np.log(empiric_mean / np.sqrt(common_val))
        sigma = np.sqrt(np.log(common_val))

        return mu, sigma

    TF = np.logspace(0, 5)
    confidence_intervals_width = [0.25, 0.5, 0.75]

    mean_vals_P = np.empty((len(TF)))
    confidence_intervals_P = np.empty((len(TF), len(confidence_intervals_width), 2))
    for iX, tf in enumerate(TF):
        m, v = model.mean_P(tf), model.var_P(tf)
        mu, sigma = log_normal_params(empiric_mean=m, empiric_var=v)
        cur_log_norm = sts.lognorm(s=sigma, scale=np.exp(mu))

        mean_vals_P[iX] = m

        for iC, confidence in enumerate(confidence_intervals_width):
            confidence_intervals_P[iX][iC] = cur_log_norm.interval(confidence=confidence)

        pass

    colors_to_use = ["#9cd3f2", "#b4def5", "#d4ecf9"]

    plt.figure(figsize=(3.8, 2.5))
    fig = plt.gcf()
    ax = plt.gca()

    # ax.set_title("Protein")

    for iI1 in range(len(confidence_intervals_width)):
        iI = len(confidence_intervals_width) - 1 - iI1
        confidence = confidence_intervals_width[iI]
        ax.fill_between(TF, confidence_intervals_P[:, iI, 1], confidence_intervals_P[:, iI, 0],
                        color=colors_to_use[iI],  # color_palette[-1],
                        # alpha=0.25,
                        label=f"{(int(confidence * 100))} %")
    # ax.plot(C, meanP_dist)
    ax.plot(TF, mean_vals_P, "--", c="magenta", label=r"$\bar{p}$")

    # ax.set_yticks([iX * 500 for iX in range(5)])
    # ax.set_yticks([iX * 4000 for iX in range(4)])
    ax.legend(prop={'size': 12})
    ax.set_xscale("log")
    # ax.set_yscale("log")

    ax.set_xlabel("Input TF")
    ax.set_ylabel("Output TF")

    plt.tight_layout()
    plt.savefig(f"protein_distribution.png", dpi=600)
    plt.show()


def visualize_protein_distribution_condensed(ax, gate, color=None):
    if color is None:
        color = "green"

    def log_normal_params(empiric_mean, empiric_var):
        common_val = empiric_var / empiric_mean ** 2 + 1
        mu = np.log(empiric_mean / np.sqrt(common_val))
        sigma = np.sqrt(np.log(common_val))

        return mu, sigma

    TF = np.logspace(0, 4)
    confidence_intervals_width = [0.25, 0.5, 0.75]

    mean_vals_P = np.empty((len(TF)))
    confidence_intervals_P = np.empty((len(TF), len(confidence_intervals_width), 2))
    for iX, tf in enumerate(TF):
        m, v = gate.model.mean_P(tf), gate.model.var_P(tf)
        mu, sigma = log_normal_params(empiric_mean=m, empiric_var=v)
        cur_log_norm = sts.lognorm(s=sigma, scale=np.exp(mu))

        mean_vals_P[iX] = m

        for iC, confidence in enumerate(confidence_intervals_width):
            confidence_intervals_P[iX][iC] = cur_log_norm.interval(confidence=confidence)

        pass

    # colors_to_use = ["#9cd3f2", "#b4def5", "#d4ecf9"]

    # plt.figure(figsize=(3.5, 2.2))
    # fig = plt.gcf()
    # ax = plt.gca()

    # ax.set_title("Protein")

    for iI1 in range(len(confidence_intervals_width)):
        iI = len(confidence_intervals_width) - 1 - iI1
        confidence = confidence_intervals_width[iI]
        ax.fill_between(TF, confidence_intervals_P[:, iI, 1], confidence_intervals_P[:, iI, 0],
                        color=color,  # color_palette[-1],
                        alpha=0.25,
                        label=f"{(int(confidence * 100))} %")
    # ax.plot(C, meanP_dist)
    ax.plot(TF, mean_vals_P, "-", c=color, label=r"$\bar{p}$")

    # ax.set_yticks([iX * 500 for iX in range(5)])
    # ax.set_yticks([iX * 4000 for iX in range(4)])
    # ax.legend()

    limits = (10 ** 0, 0.8 * 10 ** 4)
    ax.set_xlim(limits)
    ax.set_ylim(limits)

    ax.set_xscale("log")
    # ax.set_yscale("log")

    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    # plt.tick_params(
    #     axis='x',  # changes apply to the x-axis
    #     which='both',  # both major and minor ticks are affected
    #     bottom=False,  # ticks along the bottom edge are off
    #     top=False,  # ticks along the top edge are off
    #     labelbottom=False)  # labels along the bottom edge are off

    # ax.set_xlabel("Input TF")
    # ax.set_ylabel("Output TF")

    # ax.axis("off")

    # plt.tight_layout()
    # plt.savefig(f"gate_protein_distribution_simplified.png", dpi=600)
    # plt.show()


def plot_distributions_in_gatelib(gatelib):
    gate_names = list(gatelib.gates_by_type_and_name["NOT"].keys())

    fig, axes = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True, figsize=(6.299214724409, 3.93701))

    COLORS = list(sns.color_palette("hls", 12))

    for iR in range(3):
        for iC in range(4):
            iI = iR * 4 + iC
            ax = axes[iR, iC]

            gate_name = gate_names[iI]
            gate = gatelib.gates_by_type_and_name["NOT"][gate_name]

            visualize_protein_distribution_condensed(ax, gate, color=COLORS[iI])

    plt.tight_layout()
    plt.savefig("gatelib_visualized.png", dpi=600, transparent=False)
    plt.show()


def visualize_promoter_distribution(model, tf):
    distribution = model.promoter_model.distribution(tf)
    positions = np.arange(len(distribution))

    plt.figure(figsize=(0.5, 0.5))
    ax = plt.gca()

    ax.bar(positions[:2], distribution[:2], color="#E6001A")
    ax.bar(positions[2:], distribution[2:], color="#00B050")

    # ax.set_xlabel("State")

    ax.set_ylim([0, 1])
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(f"promoter_distribution_{tf}.png", dpi=600, transparent=True)
    plt.show()

    # colors = ["#E6001A", "#00B050"]
    # for iS, prob in enumerate(distribution):
    #     plt.figure(figsize=(1, 1))
    #     ax = plt.gca()
    #
    #     ax.bar(0, prob, color=colors[0] if iS < 2 else colors[1])
    #
    #     ax.set_ylim([0, 1.1])
    #     ax.axis("off")
    #     plt.tight_layout()
    #     plt.savefig(f"promoter_distribution_state_{iS}_{tf}.png", dpi=600)
    #     # plt.show()


def visualize_gate_samples(gate, tf, N=2000):
    sim_settings = {"mode": "samp"}
    samples = [gate(tf, sim_settings) for iX in range(N)]

    plt.figure(figsize=(3, 1.2))

    ax = plt.gca()
    colors = ["#9cd3f2ff", "#f29bf2ff"]
    bins = np.linspace(3 * 10 ** 3, 7 * 10 ** 3, 40)
    ax.hist(samples, bins=bins, density=True, color=colors[1])
    # ax.set_xscale("log")
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(f"simple_protein_distribution.png", dpi=600)
    plt.show()


def get_gate_by_name(gate_lib_json, name):
    relevant_gates = [gate for gate in gate_lib_json.data if gate["identifier"] == name]
    if len(relevant_gates) == 0:
        raise Exception(f"Gate {name} not in \n{gate_lib_json.content}")

    return relevant_gates[0]


if __name__ == '__main__':
    gate_lib_path = "data/gate_libs/gate_lib_yeast_cello_data_fixed_limits.json"
    gate_lib_json = JsonFile(gate_lib_path)

    gate_lib = GateLib(gate_lib_json)
    gate = gate_lib.gates_by_type_and_name["NOT"]["P1_BM3RI"]

    # visualize_gate_samples(gate, 10 ** 2)

    plot_distributions_in_gatelib(gate_lib)

    gate = get_gate_by_name(gate_lib_json, name="P1_BM3RI")

    model_params = gate["biorep"]["model"]

    model = FourStatePromoterModel(model_params)

    visualize_protein_distribution(model)

    visualize_protein_distribution_condensed(model)

    TF = np.power(10, [1, 2, 3])

    for tf in TF:
        visualize_promoter_distribution(model, tf)

    pass
