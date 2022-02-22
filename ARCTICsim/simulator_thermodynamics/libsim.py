# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 13:52:04 2020

@author: user
"""

import numpy as np
import csv
import math
import matplotlib.pyplot as plt
import json
import time
import os
import sys
from copy import copy, deepcopy

# The interfacing in the simulator methods and marshalling and demarshalling
# of objects is based on json. Therefor, all classes can be generated from
# files or strings handed over via command line. To make the objects writable
# they must implement a __str__() method for marshalling that converts them
# to a json representation, which may then be stored in a file.

# List of content
#     - technology mapping interface related
#         - classes:
#             - simulator_settings
#             - runtime_interface
#     - circuit related classes and functions
#         - classes:
#             - nand_circuit
#             - nor_circuit
#             - circuit_structure
#             - circuit_assignment
#             - library
#         - functions:

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def hl(word):
    return bcolors.WARNING + str(word) + bcolors.ENDC

#____________________________________________________________________________
#   Technology mapping interface related classes
#____________________________________________________________________________

# class simulator_settings(dict):
#     def __init__(self, *arg, **kw):
#       super(simulator_settings, self).__init__(*arg, **kw)

class execution_wrapper:
    def __init__(self, command_dict):
        self.main = command_dict
        self.pre = dict({k: (lambda *args, **kwargs: (args, kwargs)) for k, v in self.main.items()})
        self.post = dict({k: (lambda ret: ret) for k, v in self.main.items()})
        self.compose = dict({k: [True, True] for k in self.main.keys()})
    def register_pre(self, key, key_or_fun, compose=True):
        if isinstance(key_or_fun, str) and key_or_fun in self.main.keys():
            self.pre[key] = self.main[key_or_fun]
        elif callable(key_or_fun):
            self.pre[key] = fun
        else:
            raise ValueError('Supplied argument is neither callable nor a registered function.')
        self.compose[key][0] = compose
    def register_post(self, key, key_or_fun, compose=True):
        if isinstance(key_or_fun, str) and key_or_fun in self.main.keys():
            self.post[key] = self.main[key_or_fun]
        elif callable(key_or_fun):
            self.post[key] = fun
        else:
            raise ValueError('Supplied argument is neither callable nor a registered function.')
        self.compose[key][1] = compose
    def call(self, key, *args, **kwargs):
        if key not in self.main:
            raise ValueError('Function name to be called unknown.')
        if self.compose[key][0]:
            main_args, main_kwargs = self.pre[key](*args, **kwargs)
            ret = self.main[key](*main_args, **main_kwargs)
        else:
            self.pre[key]()
            ret = self.main[key](*args, **kwargs)
        if self.compose[key][1]:
            return self.post[key](ret)
        return self.post[key]()



class communication_wrapper:
    def __init__(self, input_file, output_file, prefix=None):
        self.i = None
        self.o = None
        self.prefix = prefix
        if input_file is None:
            self.i = sys.stdin
        else:
            pass # stored and temporary files not supported yet
        if output_file is None:
            self.o = sys.stdout
        else:
            pass # stored and temporary files not supported yet
    def readline(self):
        # first signal, that you are ready, by printing the prefix
        if self.prefix is not None:
            self.o.write(self.prefix + ' ')
            self.o.flush()
        return self.i.readline()
    def writeline(self, line):
        self.o.write(line + '\n')

#____________________________________________________________________________
#   Circuit related classes and functions
#____________________________________________________________________________

# input encoding dict
_input_encoding = dict({
    'a': 0,
    'b': 1,
    'c': 2
})

# A structure for a NAND circuit
class nand_circuit:
    def __init__(self, structure, assignment, library, solver=None):
        self.structure = structure
        self.assignment = assignment
        self.lib = library
        self.solver = solver
        self.promoters = dict()
        self.gates = list()
        # get a list of the involved tfs
        self.tfs = dict()
        self.tfs['internal'] = set([lib.devices[dev]['tf'] for dev in self.assignment.map_dtog.keys() if dev[0] != '_'])
        self.tfs['input'] = set([lib.devices[dev]['tf'] for dev in self.assignment.map_dtog.keys() if dev[0] != '_'])
        self.tfs['output'] = set([lib.devices[dev]['tf'] for dev in self.assignment.map_dtog.keys() if dev[0] != '_'])
        self.counts = dict_like(tfs, fill=0.)
        # first, establish a global order over all tfs and their respective gates in the circuit
        # generate promoter objects
        for prom in self.assignment.promoters:
            e_tf = dict()
            ltfe = self.lib.parts['promoters'][prom]['e']
            e_tf['rnap'] = ltfe['rnap']
            e_tf['tf_only'] = filter_dict(ltfe['tf_only'], self.tfs)
            e_tf['tf_rnap'] = filter_dict(ltfe['tf_rnap'], self.tfs)
            e_env = dict()
            lee = self.lib.env['e']
            e_env['rnap_background'] = lee['rnap_background']
            e_env['tf_background'] = filter_dict(lee['tf_background'], self.assignment.tfs)
            self.promoters[prom] = _promoter(prom, {**e_tf, **e_env, **{'beta': self.lib.env['beta']}})
        # generate gate objects, so we can fornulate the optimization function
        for gate_name in self.structure.gates + self.structure.outputs:
            g_inputs = list()
            for edge in self.structure.edges:
                if edge['target'] == gate_name:
                    prom_name = self.assignment.map_gtod[edge['source']]
                    g_inputs += [self.promoters[prom_name]]
            self.gates += [_nand_gate(gate_name, g_inputs, self.lib.parts['tfs'][lib.devices[self.assignment.map_gtod]['tf'][gate_name]]['typical']['on'])]
    def set_solver(self, solver_class):
        self.solver = solver_class(None, self.tfs, self.gates)
    def solve(self):
        if self.solver is not None:
            self.solver.solve()
            return True
        return False
    def is_valid(self):
        pass

# A structure for a NOR circuit
class OLD_nor_circuit:
    def __init__(self, structure, assignment, library, solver=None):
        self.structure = structure
        self.assignment = assignment
        self.lib = library
        self.promoters = dict()
        self.gates = list()
        # get a list of the involved tfs
        self.tfs_flat = list()
        # first, establish a global order
        self.devices = sorted(self.assignment.map_dtog.keys(), key=str.casefold)
        self.indices = dict()
        self.indices['input'] = []
        self.indices['output'] = []
        self.indices['internal'] = []
        self.all_indices = dict()
        # sort present tf's in input, internal and output
        # this is inefficient, find a better order-preserving solution
        for nn in range(len(self.devices)):
            dev = self.devices[nn]
            if self.assignment.map_dtog[dev] in self.structure.inputs:
                self.indices['input'] += [nn]
            elif self.assignment.map_dtog[dev] in self.structure.outputs:
                self.indices['output'] += [nn]
            else:
                self.indices['internal'] += [nn]
            self.all_indices[dev] = nn
            if dev[0] != '_':
                self.tfs_flat += [self.lib.devices[dev]['tfs'][0]]
        self.values_flat = np.zeros(len(self.tfs_flat))
        self.get_energies()
        for dev in self.devices:
            # Find better solution to recognize YFP
            if not self.lib.devices[dev]['promoters']:
                self.gates += [_implicit_or_gate(dev)]
            else:
                self.gates += [_nor_gate(dev, self.energies[dev], self.energies['env'], self.levels, self.lib.env['reservoir'], self.lib.env['beta'])]
        self.w = np.zeros([len(self.devices), len(self.devices)])
        for gate_in, v in self.structure.adjacency['in'].items():
            i = self.assignment.map_gtod[gate_in]
            for gate_out in v:
                o = self.assignment.map_gtod[gate_out]
                self.w[self.all_indices[i], self.all_indices[o]] = 1
        self.solver = solver(self.values_flat[self.indices['input']], self.indices['input'], self.values_flat, self.gates, self.w)
        self.substitution_mode_on = False
    def get_energies(self):
        self.energies = dict()
        self.levels = np.zeros([len(self.devices), 2])
        for nn in range(len(self.devices)):
            dev = self.devices[nn]
            # is usually just one promoter per device
            self.levels[nn] = np.array([self.lib.parts['tfs'][self.lib.devices[dev]['tfs'][0]]['typical']])
            self.energies[dev] = dict()
            self.energies[dev]['p'] = list()
            self.energies[dev]['j'] = list()
            self.energies[dev]['f'] = list()
            for prom in self.lib.devices[dev]['promoters']:
                ltfe = self.lib.parts['promoters'][prom]['e']
                self.energies[dev]['p'] = np.array([ltfe['rnap']])
                self.energies[dev]['f'] = np.array(_filter_dict(ltfe['tf_only'], self.tfs_flat, as_list=True))
                self.energies[dev]['j'] = np.array(_filter_dict(ltfe['tf_rnap'], self.tfs_flat, as_list=True))
        lee = self.lib.env['e']
        self.energies['env'] = dict()
        self.energies['env']['p'] = lee['rnap_background']
        self.energies['env']['f'] = np.array(_filter_dict(lee['tf_background'], self.tfs_flat, as_list=True))
    def set_initial_value(self, values, truthtable_idx=0):
        for n in range(len(self.indices['input'])):
            self.values_flat[self.indices['input'][n]] = values[n]
        if self.solver is not None:
            self.solver.update_inputs(self.values_flat[self.indices['input']])
        if self.substitution_mode_on:
            for gate in self.gates:
                gate.set_env(self.circuit.gate_truthtables[self.assignment.map_dtog[gate.name]][truthtable_idx])
    def set_solver(self, solver_class):
        self.solver = solver_class(None, self.values_flat, self.gates)
    def solve(self, tol, max_iter):
        if self.solver is not None:
            (err, iter) = self.solver.solve(tol, max_iter)
            self.values_flat = self.solver.values
            return (err, iter)
        return None
    def is_valid(self):
        pass
    def substitution_mode(sub_connections, sub_environments):
        self.substitution_mode_on = True
        sub_w = np.copy(self.w)
        # cut all wires where there shall be a substituted connection
        for k, v in sub_connections.items():
            i = self.assignment.map_gtod[k]
            o = self.assignment.map_gtod[v]
            sub_w[self.all_indices[i], self.all_indices[o]] = 0
        self.solver.substitution_mode(sub_w)
    def __str__(self):
        s = "Factors: " + str(self.tfs_flat)
        s += "\nValues: " + str(self.values_flat)
        s += "\nWiring:\n" + str(self.w)
        s += "\nGates: " + str([gate.name for gate in self.gates])
        for g in self.gates:
            s += '\n'
            s += 'Gate: ' + str(g)
        return s


# A structure for a NOR circuit
# we try it here with a fixed TF-order and all information
class nor_circuit:
    def __init__(self, structure, assignment, library, solver=None):
        self.structure = structure
        self.assignment = assignment
        self.lib = library
        self.gates = None

        # create a total order over all TF's and promoters
        # p_idx creates a dict of indices associated with the specific promoters
        # tf_idx creates a dict of index tuples associated with the specific tf's
        self.p_idx = dict(zip(self.lib.parts['promoters'].keys(), list(range(len(self.lib.parts['promoters'])))))
        self.tf_idx = dict({k: None for k in self.lib.parts['tfs'].keys()})
        self.dev_idx = dict({k: None for k in self.lib.devices.keys()})
        for dev, parts in self.lib.devices.items():
            for tf in parts['tfs']:
                self.tf_idx[tf] = tuple(self.p_idx[p] for p in parts['promoters'])
            self.dev_idx[dev] = tuple(self.p_idx[p] for p in parts['promoters'])

        # get factors for all possible promoters
        self.factors = np.zeros([2, len(self.p_idx), len(self.p_idx)])
        self.extremes = np.zeros([len(self.p_idx), 2])
        for p, d in self.lib.parts['promoters'].items():
            self.extremes[self.p_idx[p], 1] = float(d['typical']['on'])
            self.extremes[self.p_idx[p], 0] = float(d['typical']['off'])
            for tf in self.tf_idx.keys():
                self.factors[0, self.p_idx[p], self.tf_idx[tf]] = float(d['f']['rnap'])*float(d['f']['tf_rnap'][tf])
                self.factors[1, self.p_idx[p], self.tf_idx[tf]] = float(d['f']['tf_only'][tf])

        # create now a total order of all locally available gates
        node_list = list(self.structure.nodes)
        self.node_idx = dict(zip(node_list, list(range(len(node_list)))))

        # create an array which indicates which gate outputs are mutable
        self.mutable = np.ones(len(node_list))
        for gate in self.structure.inputs:
            self.mutable[self.node_idx[gate]] = 0

        # execute the assignment
        # get the maps g_p, the inverse p_g and w == g_g resolving the
        # circuit transformed mapping p_p == g_p ( w ( p_g ) )
        self.g_p = np.zeros(len(node_list), dtype=int)
        self.p_g = -np.ones(len(self.p_idx), dtype=int)
        self.dummy_idx = np.zeros(len(self.assignment.dummys), dtype=int)
        self.gates = np.array([None for _ in range(len(node_list))])  # array of objects
        for k, v in self.assignment.map_gtod.items():
            if v[0] == '_':  # dummy gate
                self.g_p[self.node_idx[k]] = -1
                self.dummy_idx[int(v[1])] = self.node_idx[k]
                # translate the dict describing the dummys behaviour into index notation
                # gate_idx -> env -> a_out, p_idx
                codebook = np.zeros([len(node_list), 2, 2])
                for l, w in self.assignment.dummys[v].items():
                    # TODO: this can be shorter
                    nidx = self.node_idx[l]
                    codebook[nidx, 0, 0] = float(w['a'][0])
                    codebook[nidx, 1, 0] = float(w['a'][1])
                    codebook[nidx, 0, 1] = self.p_idx[w['b'][0]]
                    codebook[nidx, 1, 1] = self.p_idx[w['b'][1]]
                self.gates[self.node_idx[k]] = _dummy_gate(k, v, codebook)
            elif v == 'OR_IMPL' or len(self.dev_idx[v]) == 0:  # implicit or
                self.g_p[self.node_idx[k]] = -1
                self.gates[self.node_idx[k]] = _implicit_or_gate(k, v)
            else:
                self.g_p[self.node_idx[k]] = self.dev_idx[v][0]
                self.p_g[self.dev_idx[v][0]] = self.node_idx[k]
                # Find better solution to recognize YFP
                self.gates[self.node_idx[k]] = _nor_gate(k, v, self.factors[:, self.dev_idx[v][0], :], self.extremes[self.dev_idx[v][0], :], self.lib.env['reservoir'])
        # also always create the artificial output gates (which are not in assignment)
        for k in list(self.structure.outputs):
            self.g_p[self.node_idx[k]] = -1
            self.gates[self.node_idx[k]] = _implicit_or_gate(k, 'simulation_monitor')

        # generate the adjacency matrix w == g_g
        self.w = np.zeros([len(node_list), len(node_list)])
        for gate_in, v in self.structure.adjacency['in'].items():
            f = self.node_idx[gate_in]
            for gate_out in v:
                t = self.node_idx[gate_out]
                self.w[f, t] = 1

        # for bounding, get a vector of extreme values
        self.bound_a_min = np.array([gate.min for gate in self.gates])
        self.bound_a_max = np.array([gate.max for gate in self.gates])
        self.bound_env = np.zeros(len(self.gates), dtype=int)

        # Now finally set the solver
        if solver is not None:
            self.set_solver(solver)

    def propagate(self, a):
        new_a = np.copy(a)
        wa = np.dot(self.w, a)
        p_a = np.zeros(len(self.p_idx))
        #print([(gate.node, self.node_idx[gate.node]) for gate in self.gates])
        for n in range(len(self.gates)):
            if self.g_p[n] != -1:
                p_a[self.g_p[n]] += wa[n]
        for n in range(len(self.gates)):
            if not self.mutable[n]:
                continue
            gate = self.gates[n]
            new_a[n] = self.gates[n].out(p_a)
        return new_a

    # this will be sloooooow.....
    def propagate_bound(self, a):
        new_a = np.copy(a)
        for n in range(len(self.gates)):
            if not self.mutable[n] or self.gates[n].type == -1:
                continue
            p_a = np.zeros(len(self.p_idx))
            # get extreme inputs
            a_x = np.zeros(len(self.gates))
            if self.bound_env[n] == 1:  # everything on minimum!
                a_x = np.copy(self.bound_a_min)
            else:  # everything on maximum!
                a_x = np.copy(self.bound_a_max)
            for idx in self.dummy_idx:  # set the dummy gates optimally
                a_x[idx] = self.gates[idx].out(n, self.bound_env[n])
            # do the wiring
            wa = np.dot(self.w.T, a_x)
            # but unset the own wire because we reset it later
            wa[n] = 0
            # now associate the correct TF's
            for m in range(len(self.gates)):
                if self.gates[m].type == -1:  # dummy: choose weakest or strongest TF
                    pa_idx = self.gates[m].promoter(n, self.bound_env[n])
                    p_a[pa_idx] += wa[m]
                else:
                    p_a[self.g_p[m]] += wa[m]
            # set the correct wire input
            for m in range(len(self.gates)):
                if self.w[m, n] == 1:  # we have a an input, so set
                    if self.gates[m].type == -1:  # dummy
                        wa[n] += self.gates[m].out(n, self.bound_env[n])
                    else:
                        wa[n] += a[m]
            # propagate the artificial environment through the gate
            new_a[n] = self.gates[n].out(p_a)
        return new_a

    def set_dummy_mode(self, ttix=0):
        # this mode must be set at least once if dummy gates are present
        for n in range(len(self.gates)):
            if self.gates[n].type != -1:  # no dummy
                self.bound_env[n] = self.structure.gate_truthtables[self.gates[n].node][ttix]

    def set_initial_value(self, values, ttix=-1):
        val = np.zeros(len(self.node_idx))
        for input in self.structure.inputs:
            val[self.node_idx[input]] = values[_input_encoding[input]]
        self.solver.set_initial_value(val)
        if ttix != -1:
            self.set_dummy_mode(ttix)

    def set_solver(self, solver_class):
        self.solver = solver_class(self, None, False)

    def solve(self, tol, max_iter):
        if self.solver is not None:
            (values, err, iter) = self.solver.solve(tol, max_iter)
            # currently supports only one output
            return (values[self.node_idx[list(self.structure.outputs)[0]]], err, iter)
        return None

    def is_valid(self):
        pass

    def __len__(self):
        return len(self.node_idx)

    def __str__(self):
        s = "\nWiring:\n" + str(self.w)
        s += "\nGates: " + str([gate.name for gate in self.gates])
        for g in self.gates:
            s += '\n'
            s += 'Gate: ' + str(g)
        return s


class circuit_structure:
    def __init__(self, str):
        circ = json.loads(str)
        if circ is not None:
            self.order = None
            self.truthtable = np.array(list(map(int, circ['truthtable'])))
            self.inputs = set()
            self.gates = set()
            self.outputs = set()
            self.nodes = set()
            self.stats = {'n_inputs': 0, 'n_nodes': 0, 'n_outputs': 0}
            for n in circ['graph']['nodes']:
                if n['type'] == 'INPUT':
                    self.stats['n_inputs'] += 1
                    self.inputs.add(n['id'])
                elif n['type'] == 'OUTPUT':
                    self.stats['n_outputs'] += 1
                    self.outputs.add(n['id'])
                else:
                    self.stats['n_nodes'] += 1
                    self.gates.add(n['id'])
            self.nodes = self.inputs | self.gates | self.outputs
            self.gate_truthtables = dict()
            for k, v in circ['gate_truthtables'].items():
                self.gate_truthtables[k] = np.array(list(map(int, v)))
            self.internal_edges = set()
            self.outgoing_edges = set()
            self.adjacency = {'in': dict_from_list(self.nodes, set()), 'out': dict_from_list(self.nodes, set())}
            for e in circ['graph']['edges']:
                if e['source'] in self.gates and e['target'] in self.gates:
                    self.internal_edges.add(tuple((e['source'], e['target'])))
                else:
                    self.outgoing_edges.add(tuple((e['source'], e['target'])))
                self.adjacency['in'][e['target']].add(e['source'])
                self.adjacency['out'][e['source']].add(e['target'])
                # if e['source'] in self.inputs:
                #     self.adjacency['in'][e['target']].add(e['source'])
                # elif e['target'] in self.outputs:
                #     self.adjacency['out'][e['source']].add(e['target'])
            self.edges = set.union(self.internal_edges, self.outgoing_edges)
            #self.dprint('\n\t+ Circuit has truthtable {f}, {u} inputs, {y} outputs, {n} intermediate nodes as well as {e} internal and {o} outgoing edges.'.format(f=hl(self.truthtable), u=hl(self.stats['n_inputs']), y=hl(self.stats['n_outputs']), n=hl(self.stats['n_nodes']), e=hl(len(self.internal_edges)), o=hl(len(self.outgoing_edges))))
            self.valid = True
    def combinational_order(self):
        if self.order is None:
            self.order = list()
            for gate in self.gates:
                for n in range(len(self.order)):
                    if self.order[n] in self.adjacency['out'][gate]:
                        self.order.insert(n, gate)
        return self.order
    # return the JSON representation of this graph
    def __str__(self):
        pass

# the assignment connects the library with the structure. It therefore needs to know both to sanity check
class circuit_assignment:
    def __init__(self, jstr):
        assi = json.loads(jstr)
        self.map_gtod = dict()
        self.map_dtog = dict()
        self.dummys = dict()
        n_dummy = 0
        for k, v in assi.items():
            d = v
            if isinstance(v, dict):
                # if there is a dict present, we want a dummy gate which emits extreme values
                d = '_' + str(n_dummy)
                self.dummys[d] = v
                n_dummy += 1
            self.map_gtod[k] = d
            self.map_dtog[d] = k
        self.devices = set(self.map_dtog.keys())
        self.gates = set(self.map_dtog.values())

class library:
    # represents a connected gate lib
    def __init__(self, jstr):
        self.parts = dict()
        self.parts['promoters'] = dict()
        self.parts['tfs'] = dict()
        self.devices = dict()
        self.env = dict()
        libstr = json.loads(jstr)
        for entry in libstr:
            if entry['class'] == 'transcription_factors':
                for part in entry['members']:
                    self.parts['tfs'][part['name']] = dict()
                    self.parts['tfs'][part['name']]['part_of'] = part['associated_devices']
                    for dev in part['associated_devices']:
                        if dev not in self.devices:
                            self.devices[dev] = dict()
                            self.devices[dev]['tfs'] = list()
                            self.devices[dev]['promoters'] = list()
                        self.devices[dev]['tfs'] += [part['name']]
                    self.parts['tfs'][part['name']]['typical'] = [part['levels']['off'], part['levels']['on']]
            if entry['class'] == 'promoters':
                for part in entry['members']:
                    self.parts['promoters'][part['name']] = dict()
                    self.parts['promoters'][part['name']]['part_of'] = part['associated_devices']
                    for dev in part['associated_devices']:
                        self.devices[dev]['promoters'] += [part['name']]
                    if 'levels' in part:
                        self.parts['promoters'][part['name']]['typical'] = part['levels']
                    self.parts['promoters'][part['name']]['e'] = part['energies']
                    self.parts['promoters'][part['name']]['f'] = part['factors']
            if entry['class'] == 'environment':
                self.env['name'] = entry['name']
                self.env['beta'] = entry['beta']
                self.env['reservoir'] = entry['non_specific_reservoir']
                self.env['df'] = entry['env_tfs']
                self.env['e'] = entry['reservoir_energies']
    # return the JSON representation of this lib
    def __str__(self):
        content = list()
        content += [dict()]
        content[0]['class'] = 'transcription_factors'
        content[0]['members'] = list()
        for k, v in self.parts['tfs'].items():
            entry = dict()
            entry['name'] = k
            entry['associated_devices'] = v['part_of']
            entry['levels'] = dict()
            entry['levels']['off'] = v['typical'][0]
            entry['levels']['on'] = v['typical'][1]
            content[0]['members'] += [entry]
        content += [dict()]
        content[1]['class'] = 'promoters'
        content[1]['members'] = list()
        for k, v in self.parts['promoters'].items():
            entry = dict()
            entry['name'] = k
            entry['associated_devices'] = v['part_of']
            entry['energies'] = v['energies']
            content[1]['members'] += [entry]
        content += [dict()]
        content[2]['class'] = 'environment'
        content[2]['name'] = self.env['name']
        content[2]['beta'] = self.env['beta']
        content[2]['non_specific_reservoir'] = self.env['reservoir']
        content[2]['env_tfs'] = self.env['df']
        content[2]['reservoir_energies'] = self.env['e']
        return json.dumps(content)

# PRIVATE. Not intended to use manually
class _promoter:
    # represents the promoter model
    def __init__(self, name, e):
        self.name = name
        self.ep = np.exp(-e['beta']*(e['rnap'] - e['rnap_background']))
        self.ef = np.exp(-e['beta']*(np.array(list(e['tf_only'].values())) - np.array(list(e['tf_background'].values()))))
        self.ej = np.exp(-e['beta']*(np.array(list(e['tf_rnap'].values())) - e['rnap'] - np.array(list(e['tf_background'].values()))))

# PRIVATE. Not intended to use manually
class _nand_gate:
    # represents a NAND gate in the circuit
    def __init__(self, name, promoters, b):
        self.name = name
        self.p = promoters
        self.w = []
        acc = 0
        for p in promoters:
            self.w += [p.ep]
            acc += p.ep
        self.w = np.array(self.w)/acc
    def out(self, f):
        # input g is a gate structure, which contains a set of promoters and an output TF
        return np.dot(self.w, np.array([_D(p, f) for p in self.p]))*self.b
    def _D(self, f):
        # The D are the inidividual regulation factors of the promoters
        return (1 + np.sum(p.ej*f))/(1 + np.sum(p.ef*f))

# PRIVATE. Not intended to use manually
# The NOR gate is a simple, modular gate structure independent of input wiring
class OLD_nor_gate:
    # represents a NOR gate in the circuit
    def __init__(self, name, e_local, e_env, b, c, beta):
        self.name = name
        # energies = {p = np.array(tfs ->), f = np.array(tfs ->), j = np.array(tfs ->)}
        self.e = dict()
        self.e['p'] = np.exp(-beta*(e_local['p'] - e_env['p']))
        self.e['f'] = np.exp(-beta*(e_local['f'] - e_env['f']))
        self.e['j'] = np.exp(-beta*(e_local['j'] - e_local['p'] - e_env['f']))
        # b = np.array(tfs ->)
        self.b = b[:, 1]
        # precalculate all the stuff that can be calculated
        self.bepj = self.b*self.e['p']*self.e['j']
        self.bef = self.b*self.e['f']
        self.c = c
    def out(self, wa):
        return (self.c + np.sum(wa*self.bepj))/(self.c + np.sum(wa*self.bef))
    def __str__(self):
        return self.name + ': ' + str(self.e)
    def set_env(self, env):
        self.env = env

# PRIVATE. Not intended to use manually
# The NOR gate is a simple, modular gate structure independent of input wiring
class _nor_gate:
    # represents a NOR gate in the circuit
    def __init__(self, node, dev, factors, extreme, c):
        self.node = node
        self.dev = dev
        self.min = extreme[0]
        self.max = extreme[1]
        self.bepj = factors[0]
        self.bef = factors[1]
        self.c = c
        self.type = 0
    def out(self, wa):
        return (self.c + np.sum(wa*self.bepj))/(self.c + np.sum(wa*self.bef))
    def __str__(self):
        return self.name + ': ' + str(self.e)

# PRIVATE. Not intended to use manually
# The implicit OR gate is just a pseudo-gate, which outputs the YFP RPU
class _implicit_or_gate:
    # represents a NOR gate in the circuit
    def __init__(self, node, dev):
        self.node = node
        self.dev = dev
        self.type = 1
        self.min = 0
        self.max = 0
    def out(self, wa):
        return np.sum(wa)
    def __str__(self):
        return self.name + ': None (implicit OR)'

# PRIVATE. Not intended to use manually
# The dummy gate is a placeholder, which always emits 0 but has initialised
# boundary values
class _dummy_gate:
    # represents a NOR gate in the circuit
    def __init__(self, node, dev, codebook):
        self.node = node
        self.dev = dev
        self.codebook = codebook
        self.type = -1
        self.min = 0
        self.max = 0
    def out(self, g_idx, env):
        return self.codebook[g_idx, env, 0]
    def promoter(self, g_idx, env):
        return int(self.codebook[g_idx, env, 1])
    def __str__(self):
        return self.name + ': None (dummy gate, codebook size: ' + len(self.codebook[:, 0, 0]) + ')'

# Circuit solver classes. Pass any of your choice to a circuit object
# This solver does a fixed point iteration. Simplest solver thinkable
class OLDnor_circuit_solver_banach:
    def __init__(self, inputs, input_idx, values, gates, wiring):
        self.inputs = inputs
        self.input_idx = input_idx
        self.W = np.transpose(wiring)
        self.values = values
        self.gates = gates
        self.substitution_mode_on = False
    def update_inputs(self, inputs):
        self.inputs = inputs
        print(inputs)
    def solve(self, tol=10**(-2), max_iter=100):
        # Iterate the points using Banach's fixed point theorem
        # to estimate the error, compare old and new node values
        err = np.inf
        run = 0
        if self.substitution_mode_on:
            wiring = self._substitution_wiring
        else:
            wiring = self._exact_wiring
        vs = self.values
        while err > tol:
            new_vs = np.array([g.out(wiring(vs, g)) for g in self.gates])
            new_vs[self.input_idx] = np.copy(self.inputs)
            # error and max_iter updates
            err = np.max(np.abs(new_vs - vs))
            vs = new_vs
            run += 1
            if run >= max_iter:
                break
        self.values = vs
        return (err, run)
    def substitution_mode(sub_wiring):
        # get the extreme values of all gates
        self.substitution_mode_on = True
        self.sub_values = np.zeros([len(self.values), 2])
        nn = 0
        for g in self.gates:
            self.sub_values[nn, :] = g.minmax
            nn += 1
        self.sub_W = np.transpose(sub_wiring)
    def _exact_wiring(self, vs, g):
        return np.dot(self.W, vs)
    def _substitution_wiring(self, vs, g):
        # exchange all values with most beneficial possible (guarantee environment)
        r_vs = np.copy(self.sub_values[g.sub_env])
        # realize circuit and guarantee valid upper bound (only 1-/0- inputs)
        r_vs[self.sub_W[n]] = np.array([self.__sub_extreme(vs[self.sub_W[n]], g.sub_env) for n in range(len(self.sub_W))])
        return r_vs
    # static method
    def __sub_extreme(data, env):
        if env == 1:
            return np.min(data)
        return np.max(data)

# This solver does a fixed point iteration. Simplest solver thinkable
# this dramatically reduced. Essentially only solve the propagate function
class nor_circuit_solver_banach:
    def __init__(self, circuit, values, bound=False):
        self.circuit = circuit
        if values is not None:
            self.values = np.copy(values)
        else:
            self.values = np.zeros(len(self.circuit))
        self.bound = bound
    def set_initial_value(self, initial_value):
        self.values = initial_value
    def bounding_mode(self, bound):
        self.bound = bound
    def solve(self, tol=10**(-2), max_iter=100):
        # Iterate the points using Banach's fixed point theorem
        # to estimate the error, compare old and new node values
        err = np.inf
        run = 0
        if self.bound:
            function = self.circuit.propagate_bound
        else:
            function = self.circuit.propagate
        vs = self.values
        while err > tol:
            new_vs = function(vs)
            # error and max_iter updates
            err = np.max(np.abs(new_vs - vs))
            vs = new_vs
            run += 1
            if run >= max_iter:
                break
        self.values = vs
        return (vs, err, run)

# Circuit solver classes. Pass any of your choice to a circuit object
# This solver does a fixed point iteration. Simplest solver thinkable
class nand_circuit_solver_banach:
    def __init__(self, inputs, tfs, gates):
        self.inputs = inputs
        self.tfs = tfs
        self.gates = gates
    def solve(self, tol=10**(-2), max_iter=100):
        # Iterate the points using Banach's fixed point theorem
        # to estimate the error, compare old and new node values
        err = np.inf
        run = 0
        while err > tol:
            new_tfs = [g.out(tfs) for g in self.gates]
            # error and max_iter updates
            err = np.max(np.abs(new_tfs - tfs))
            tfs = new_tfs
            run += 1
            if run >= max_iter:
                break
        return (err, run)

# This solver uses Newton's method. Reliable and fast (quadratic convergence)
class circuit_solver_newton:
    def __init__(self, inputs, tfs, gates):
        pass

#____________________________________________________________________________
#   File interfacing classes
#____________________________________________________________________________

class json_file:
    def __init__(self, path, read=True):
        self.path = path
        self.content = None
        if read:
            self.read()
    def read(self, update=False):
        if self.content is None or update == True:
            with open(self.path, 'r') as fyle:
                self.content = fyle.read().replace('\n', '')#json.load(fyle)
    def write(self):
        with open(self.path, 'w') as fyle:
            self.content = fyle.write(self.content)#json.dump(fyle)

#____________________________________________________________________________
#   Classes and functions for more abstract objects/tasks
#____________________________________________________________________________

class __directed_acyclic_graph:
    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges
    def is_valid(self):
        pass

def __list_intersect(a, b):
    return list(set(a) & set(b))

def _filter_dict(d, l, as_list=False, act=copy):
    if as_list:
        return [act(d[ll]) for ll in l if ll in d]
    return {ll: act(d[ll]) for ll in l if ll in d}

def dict_like(d, fill=None):
    if not isinstance(v, dict):
        return copy(fill)
    return dict({k: dict_like(v, fill) for k, v in d.items()})

def dict_from_list(l, fill=None):
    return dict({k: copy(fill) for k in l})
