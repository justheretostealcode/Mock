#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:14:26 2020

@author: nicolai
"""

import math
import numpy as np
from configparser import ConfigParser
from argparse import ArgumentParser
import os.path as op
import os
import inspect
import itertools
import shlex
import libsim_erich as sim

# standard configuration path
version = "1.2"
settings_config_path = "./settings_config_erich.cfg"

# pwd
here = op.dirname(op.abspath(inspect.getfile(inspect.currentframe())))

# ____________________________________________________________________________
#   Global variables to share among threads written only beforehand
# ____________________________________________________________________________

library = None
circuit = None
structure = None
inputs = None
supplement = None

library_file = None
structure_file = None
assignment_file = None
supplement_file = None

# ____________________________________________________________________________
#   Functions and related constants for convenience
# ____________________________________________________________________________


# load settings from a file. Descriptions of settings are also loaded
def load_settings(path):
    s = ConfigParser(allow_no_value=True)
    s.read(path)
    s = dict({k: v for sec in s.sections() for k, v in s[sec].items()})
    c = dict()
    t = dict()
    for k, v in s.items():
        s[k], t[k], c[k] = tuple(map(str.strip, v.split("#", 2)))
        s[k] = type_dict[t[k]](s[k])
    return s, t, c


# an associated type dict (for safety we do not use eval)
type_dict = dict(
    {
        "int": int,
        "str": str,
        "bool": bool,
        "float": float,
    }
)


def check_load_and_create(input, file, cluss):
    # check if what we got is a file or a marshalled json object
    if op.exists(input):
        file = sim.json_file(input)
        json = file.content
    else:
        json = input
    # demarshall the object
    return cluss(json)


def real_from_logic_input_matrix(linputs, library, circuit):
    n_ins = len(linputs[0])
    on_off = ["off", "on"]
    real_inputs = [[0.0034, 2.8], [0.0013, 4.4], [0.0082, 2.5]]
    inputs = np.zeros([len(linputs), n_ins])
    for n in range(len(linputs)):
        for m in range(n_ins):
            lin_to_rin = circuit.gates[circuit.node_idx[sim.input_encoding[m]]].dev
            if not lin_to_rin.startswith("_"):
                inputs[n][m] = real_inputs[m][int(not bool(linputs[n][m]))]
            else:  # dummy
                inputs[n][m] = 0.0
    return inputs


def round_to_magnitude(r, b):
    f = (b % 1) * 10
    m = 1
    while int(f % 10) == 0:
        m += 1
        f *= 10
    return round(r, m)


# ____________________________________________________________________________
#   Functions related to executing the actual simulation
# ____________________________________________________________________________


# check if we have all data for simulation and do some preparations
def prepare_simulation():
    global settings, library, circuit, structure, inputs, library_file, structure_file, assignment_file
    # load library
    if library is None:
        library_file = sim.json_file(settings["library"])
        library = sim.library(library_file.content)
    # check and load structure
    if structure is None:
        structure = check_load_and_create(
            settings["structure"], structure_file, sim.circuit_structure
        )
    # check and load assignment
    assignment = check_load_and_create(
        settings["assignment"], assignment_file, sim.circuit_assignment
    )
    # check and load bb additional library
    # supplement = check_load_and_create(settings['bb_library'], bb_file, sim.circuit_supplement)
    # build the circuit
    circuit = sim.nor_circuit(
        structure, assignment, library, solver=sim.nor_circuit_solver_banach
    )
    # generate the inputs (I assume, I have a valid truthtable given)
    # we further assume, the inputs are numerated in an alphabethical order
    n_ins = int(round(math.log(len(circuit.structure.truthtable), 2)))
    linputs = np.array(list(map(list, itertools.product([0, 1], repeat=n_ins))))
    inputs = real_from_logic_input_matrix(linputs, library, circuit)
    # Done with preparation. Ready for simulation now


# ____________________________________________________________________________
#   Functions and global constants related to the CLI (start, exit, update, ...)
# ____________________________________________________________________________


# Exit the simulator
def sim_exit():
    # everything here for a controlled exit
    exit(0)


# This is somehow dirty. Change that
def sim_update_settings(args):
    global settings, argp
    settings.update(
        {
            k: v
            for k, v in vars(argp.parse_args(shlex.split(args, posix=False))).items()
            if v is not None
        }
    )


# Print the current simulation settings. Again a little dirty
def sim_print_settings():
    global settings
    print(settings)


# Print the current simulation settings. Again a little dirty
def sim_print_commands():
    global command_dict
    for k, v in command_dict.items():
        print("{k:<30} -> {v}".format(k=k, v=v[1]))


# Start the simulator
def sim_run():
    global settings, circuit, inputs, assignment
    print("Running simulation.")
    # Do the preparation:
    prepare_simulation()
    bounding_mode = False
    whitelist = np.ones(len(inputs), dtype=int)
    if (
        len(circuit.assignment.dummys) > 0
        or "propagation_mode" in settings
        and settings["propagation_mode"] > 0
    ):
        bounding_mode = True
        sub_mode = settings["propagation_mode"] - 1
        circuit.solver.bounding_mode(sub_mode)
        if "whitelist" in settings:
            whitelist = np.array(list(map(int, settings["whitelist"])))
    # Now solve the circuit function
    N = np.sum(whitelist)
    N_on = np.sum(whitelist & circuit.structure.truthtable)
    results = [np.zeros(N - N_on), np.zeros(N_on)]
    used_inputs = [
        (inputs[n], circuit.structure.truthtable[n], n)
        for n in range(len(inputs))
        if whitelist[n] == 1
    ]
    c = np.zeros(2, dtype=int)
    max_err = 0.0
    for real_input, bool_output, ttix in used_inputs:
        circuit.set_initial_value(real_input, ttix)
        output, err, iter = circuit.solve(
            tol=float(settings["err"]), max_iter=int(settings["max_iter"])
        )
        results[bool_output][c[bool_output]] = round_to_magnitude(
            output, 2 * settings["err"]
        )
        max_err = max(max_err, err)
        c[bool_output] += 1
    print("results: " + str(results))
    return circuit_score(results, max_err)


# Just for now the only circuit score
# This will anyway be exchanged with a proper post-processing later
def circuit_score(results, err):
    global settings
    differences = results[1][:, None] / results[0]
    score = round(np.min(np.min(differences)), int(settings["final_precision"]))
    if err > 2 * settings["err"]:
        cli_io.writeline("error above tolerance, score dismissed.")
        return -1
    cli_io.writeline("score: " + str(score))
    return score


# to add a command to the simulator, extend this list and add your function
# any function *can* accept a string as an argument containing the further
# part of the command line
command_dict = dict(
    {
        "start": [sim_run, "starts the simulation"],
        "exit": [sim_exit, "terminates the simulator"],
        "update_settings": [sim_update_settings, "update the simulation settings"],
        "print_settings": [sim_print_settings, "print the current simulation settings"],
        "help": [sim_print_commands, "prints this help"],
    }
)

# ____________________________________________________________________________
#   The main script starts here
# ____________________________________________________________________________


if __name__ == "__main__":
    # change working directory to here
    os.chdir(here)

    # load default settings first
    settings, types, comments = load_settings(settings_config_path)

    # update settings through the command line
    argp = ArgumentParser(description="Thermodynamic circuit simulator v." + version)
    for k, v in settings.items():
        argp.add_argument(
            "-" + k[0],
            "--" + k,
            required=False,
            type=type_dict[types[k]],
            help=comments[k],
        )
    argp.add_argument(
        "--autostart",
        required=False,
        action="store_true",
        default=False,
        help="if set, simulation starts immediately (no CLI beforehand)",
    )
    argp.add_argument(
        "--autoexit",
        required=False,
        action="store_true",
        default=False,
        help="if set, termination immediately after first simulation (no CLI afterwards)",
    )
    argp.add_argument(
        "--no_cli",
        required=False,
        action="store_true",
        default=False,
        help="if set, no CLI is called at all (immediate start, immediate exit)",
    )
    settings.update({k: v for k, v in vars(argp.parse_args()).items() if v is not None})

    # setup communication interface
    cli_io = sim.communication_wrapper(None, None, prefix=":>")
    cli_exec = sim.execution_wrapper({k: v[0] for k, v in command_dict.items()})

    # set verbosity level
    sim.DEBUG_LEVEL = settings["verbosity"]

    # check for autoexit
    if settings["autoexit"] or settings["no_cli"]:
        cli_exec.register_post("start", "exit", compose=False)

    # check for autostart
    if settings["autostart"] or settings["no_cli"]:
        cli_exec.call("start")

    cli_io.writeline("ready")

    # run CLI loop
    while True:
        # read in the next line and split it into command and arguments
        line_cmd, *line_args = tuple(map(str.strip, cli_io.readline().split(" ", 1)))
        cli_exec.call(line_cmd, *line_args)
