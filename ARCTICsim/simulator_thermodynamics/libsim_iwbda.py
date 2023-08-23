# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 13:52:04 2020

@author: user
"""

import numpy as np
import numpy.linalg as nla
import csv
import math
import matplotlib.pyplot as plt
import json
import time
import os
import sys
import scipy.optimize as so
from copy import copy, deepcopy
from autograd import elementwise_grad as grad

# The interfacing in the simulator methods and marshalling and demarshalling
# of objects is based on json. Therefor, all objects can be generated from
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

#____________________________________________________________________________
#   Globals and general purpose short-cuts
#____________________________________________________________________________


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

def head(word):
    return bcolors.OKCYAN + str(word) + bcolors.ENDC

DEBUG_LEVEL = 2

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
input_encoding = dict({
    0: 'a',
    1: 'b',
    2: 'c'
})

# reverse input encoding dict
reverse_input_encoding = dict({
    'a': 0,
    'b': 1,
    'c': 2
})

# bounding submodes
bounding_modes = dict({
    'heuristic': 1,
    'ita': 2
})

# Define important indexing variables
# leaves of the decision tree
_wired_assignedandnor = 0
_wired_assignedandor = 1
_wired_dummy_nor = 2
_wired_dummy_or = 3
_xtalk_assignedandnor = 4
_xtalk_dummyoror_dummy = 5
_xtalk_dummyoror_or = 6
_xtalk_dummyoror_tf = 7
# Variable categories from the LUT's of the gates
_a = 0
_i = 1
_j = 2
_k = 3
_d = 4
_b = 5
_o = 6
_c = 7
_e = 8
_var_map = dict({
    'a': _a,
    'i': _i,
    'j': _j,
    'k': _k,
    'd': _d,
    'b': _b,
    'o': _o,
    'c': _c,
    'e': _e,
})
_var_pre = dict({
    'a': (lambda a, circ: float(a)),
    'i': (lambda i, circ: float(i)),
    'j': (lambda j, circ: float(j)),
    'k': (lambda k, circ: float(k)),
    'd': (lambda d, circ: circ.p_idx[d]),
    'b': (lambda b, circ: circ.tf_idx[b][0]),
    'o': (lambda o, circ: float(o)),
    'c': (lambda c, circ: float(c)),
    'e': (lambda e, circ: float(e)),
})
# Indexes for min and max values of a respective category
_min = 0
_max = 1

# bounding mode settings
# this defines, in which modes which variables are drawn given the parameters
_bounding_configs = dict({
    0: np.array([  # TODO: naive optimal mode
        ]),
    1: np.array([  # TODO: naive herusitic mode
        ]),
    2: np.array([  # ita optimal mode
            [  # _wired_assignedandnor
                [(_i, _min), (_o, _max)],   # 0 target
                [(_o, _min), (_i, _max)],   # 1
                #    0              1
                #       source
            ],
            [  # _wired_assignedandor
                [(_o, _min), (_i, _max)],   # 0 target
                [(_i, _min), (_o, _max)],   # 1
                #    0              1
                #       source
            ],
            [  # _wired_dummy_nor
                [(_i, _min), (_a, _max)],   # 0 target
                [(_a, _min), (_i, _max)],   # 1
                #    0              1
                #       source
            ],
            [  # _wired_dummy_or
                [(_a, _min), (_i, _max)],   # 0 target
                [(_i, _min), (_a, _max)],   # 1
                #    0              1
                #       source
            ],
            [  # _xtalk_assignedandnor
                [(_c, _max), (_j, _min)],   # 0 target
                [(_j, _max), (_c, _min)],   # 1
                #    0              1
                #       source
            ],
            [  # _xtalk_dummyoror_dummy
                [(_c, _max), (_j, _min)],   # 0 target
                [(_j, _max), (_c, _min)],   # 1
                #    0              1
                #       source
            ],
            [  # _xtalk_dummyoror_or
                [(_k, _min), (_c, _max)],   # 0 target
                [(_c, _min), (_k, _max)],   # 1
                #    0              1
                #       source
            ],
            [  # _xtalk_dummyoror_tf
                [(_b, _max), (_b, _max)],   # 0 target
                [(_b, _min), (_b, _min)],   # 1
                #    0              1
                #       source
            ],
        ]),
    3: np.array([  # TODO: ita heuristic mode
        ]),
    4: np.array([  # TODO: full heuristic mode
        ]),
    5: np.array([  # optimal mode
            [  # _wired_assignedandnor
                [(_a, _max), (_o, _max)],   # 0 target
                [(_o, _min), (_o, _max)],   # 1
                #    0              1
                #       source
            ],
            [  # _wired_assignedandor
                [(_o, _min), (_o, _max)],   # 0 target
                [(_a, _max), (_o, _max)],   # 1
                #    0              1
                #       source
            ],
            [  # _wired_dummy_nor
                [(_a, _max), (_a, _max)],   # 0 target
                [(_a, _min), (_a, _max)],   # 1
                #    0              1
                #       source
            ],
            [  # _wired_dummy_or
                [(_a, _min), (_a, _max)],   # 0 target
                [(_a, _max), (_a, _max)],   # 1
                #    0              1
                #       source
            ],
            [  # _xtalk_assignedandnor
                [(_c, _max), (_c, _max)],   # 0 target
                [(_c, _min), (_c, _min)],   # 1
                #    0              1
                #       source
            ],
            [  # _xtalk_dummyoror_dummy
                [(_c, _max), (_c, _max)],   # 0 target
                [(_c, _min), (_c, _min)],   # 1
                #    0              1
                #       source
            ],
            [  # _xtalk_dummyoror_or
                [(_c, _max), (_c, _max)],   # 0 target
                [(_c, _min), (_c, _min)],   # 1
                #    0              1
                #       source
            ],
            [  # _xtalk_dummyoror_tf
                [(_b, _max), (_b, _max)],   # 0 target
                [(_b, _min), (_b, _min)],   # 1
                #    0              1
                #       source
            ],
        ]),
})
# force maps
_force_configs = dict({
    0: np.array([  # TODO: naive optimal mode
        ]),
    1: np.array([  # TODO: naive herusitic mode
        ]),
    2: np.array([  # ita optimal mode
            [  # _wired_assignedandnor
                [(0, 0), (0, 1)],   # 0 target
                [(1, 0), (0, 0)],   # 1
                #   0       1
                #    source
            ],
            [  # _wired_assignedandor
                [(0, 0), (1, 0)],   # 0 target
                [(1, 0), (1, 1)],   # 1
                #   0       1
                #    source
            ],
            [  # _wired_dummy_nor
                [(0, 0), (0, 1)],   # 0 target
                [(1, 0), (0, 0)],   # 1
                #   0       1
                #    source
            ],
            [  # _wired_dummy_or
                [(0, 0), (1, 0)],   # 0 target
                [(1, 0), (1, 1)],   # 1
                #   0       1
                #    source
            ],
        ]),
    3: np.array([  # TODO: ita heuristic mode
        ]),
    4: np.array([  # TODO: full heuristic mode
        ]),
    5: np.array([  # optimal mode
            [  # _wired_assignedandnor
                [(0, 0), (0, 1)],   # 0 target
                [(1, 0), (0, 0)],   # 1
                #   0       1
                #    source
            ],
            [  # _wired_assignedandor
                [(0, 0), (1, 0)],   # 0 target
                [(1, 0), (1, 1)],   # 1
                #   0       1
                #    source
            ],
            [  # _wired_dummy_nor
                [(0, 0), (0, 1)],   # 0 target
                [(1, 0), (0, 0)],   # 1
                #   0       1
                #    source
            ],
            [  # _wired_dummy_or
                [(0, 0), (1, 0)],   # 0 target
                [(1, 0), (1, 1)],   # 1
                #   0       1
                #    source
            ],
        ]),
})

# A structure for a NOR circuit
# we try it here with a fixed TF-order and all information
class nor_circuit:
    def __init__(self, structure, assignment, library, solver=None):
        self.structure = structure
        self.assignment = assignment
        self.lib = library
        self.gates = None
        self.wire_capacity = 3  # change that to something adaptive

        # create a total order over all TF's and promoters
        # g_idx creates a dict of indices associated with the specific gates
        # tf_idx creates a dict of index tuples associated with the specific tf's
        #self.g_idx = dict(zip([e['name'] for e in self.lib.c['gates']], list(range(len(self.lib.c['gates'])))))
        self.dev_idx = dict(zip([e['name'] for e in self.lib.c['devices']], list(range(len(self.lib.c['devices'])))))
        self.tf_idx = dict({e['name']: tuple() for e in self.lib.c['transcription_factors']})
        #self.dev_idx = dict({e['name']: None for e in self.lib.c['devices']})
        self.g_idx = dict({e['name']: tuple() for e in self.lib.c['gates']})
        for gate in self.lib.c['gates']:
            # Assert that there is only one gate per device
            #self.dev_idx[dev['name']] = tuple(int(self.g_idx[gate['name']]) for gate in self.lib.c['gates'] if dev['name'] in gate['associated_devices'])
            self.g_idx[gate['name']] = tuple(int(self.dev_idx[dev_name]) for dev_name in gate['associated_devices'])
        for tf in self.lib.c['transcription_factors']:
            self.tf_idx[tf['name']] = tuple(int(self.dev_idx[dev_name]) for dev_name in tf['associated_devices'])

        # get parameters for all possible gates
        self.parameters = list([None for _ in range(len(self.dev_idx))])
        for gate in self.lib.c['gates']:
            self.parameters[self.g_idx[gate['name']][0]] = gate['parameters']

        # create now a total order of all locally available gates
        node_list = list(self.structure.nodes)
        self.node_idx = dict(zip(node_list, list(range(len(node_list)))))

        # create an array which indicates which gate outputs are mutable
        self.mutable = np.ones(len(node_list))
        for input in self.structure.inputs:
            self.mutable[self.node_idx[input]] = 0

        # execute the assignment
        # get the maps g_p, the inverse p_g and w == g_g resolving the
        # circuit transformed mapping p_p == g_p ( w ( p_g ) )
        self.g_p = np.zeros(len(node_list), dtype=int)
        self.p_g = -np.ones(len(self.dev_idx), dtype=int)
        self.assigned = np.ones(len(node_list))
        self.gates = np.array([None for _ in range(len(node_list))])  # array of objects
        for k, v in self.assignment.map_gtod.items():
            #print('node/device: ' + hl(str(k)) + '/' + hl(str(v)) + ', node idx = ' + hl(str(self.node_idx[k])) + ', device idx = ' + hl(str(self.dev_idx[v])))
            self.g_p[self.node_idx[k]] = self.dev_idx[v]
            self.p_g[self.dev_idx[v]] = self.node_idx[k]
            if v.startswith('output'):  # implicit or / reporter TF
                self.gates[self.node_idx[k]] = _feedback_implicit_or_gate(k, v, self.parameters[self.dev_idx[v]])
            elif v.startswith('input'):  # an input
                self.gates[self.node_idx[k]] = _base_gate(k, v)
            else:  # assigned NOR gate
                self.gates[self.node_idx[k]] = _feedback_nor_gate(k, v, self.parameters[self.dev_idx[v]])
        if True:#DEBUG_LEVEL > 0:
            print('Gates comprising the circuit:')
            for gate in self.gates:
                print(hl(str(gate.node)) + '/' + hl(str(gate.dev)) + ', type = ' + str(gate.type))

        # generate the adjacency matrix w == g_g
        # also inject a list of inputs to each gate
        if DEBUG_LEVEL > 0:
            print('Drawing wires:')
        self.w = np.zeros([len(node_list), len(node_list)])
        for gate_in, v in self.structure.adjacency['in'].items():
            t = self.node_idx[gate_in]
            in_list = list()
            for gate_out in v:
                f = self.node_idx[gate_out]
                if DEBUG_LEVEL > 0:
                    print(hl(gate_out) + ' -> ' + hl(gate_in))
                self.w[f, t] = 1
                in_list.append(f)
            self.gates[t].inputs = np.array(in_list)

        self.bound_env = np.zeros(len(self.gates), dtype=int)

        # Now finally set the solver
        if solver is not None:
            self.set_solver(solver)

    # a is the array of current output values
    # TODO: modify for a general wire capacity
    # a has three entries here:
    # a[:, 0] = vector of promoter activities
    # a[:, 1] = vector of RNA concentrations
    # a[:, 2] = vector of RNA concentrations times their local delta
    def propagate(self, a):
        new_a = np.empty_like(a)
        for n in range(len(self.gates)):
            if not self.mutable[n]:
                new_a[n] = np.copy(a[n])
                continue
            gate = self.gates[n]
            # BEGIN: this block shall later become the wire function
            input = np.empty(self.wire_capacity)
            input[0] = np.dot(self.w[:, n], a[:, 0])
            # temporary: shift delta to outgoing gates
            input[1] = np.dot(self.w[:, n], a[:, 1])
            input[2] = np.dot(self.w[n, :], a[:, 2])
            # END: this block shall later become the wire function
            new_a[n] = gate.out(input)
        return new_a

    def set_truthtable(self, ttix=0):
        # this mode must be set at least once if dummy gates are present
        if self.structure.gate_truthtables is None:
            if DEBUG_LEVEL > 0:
                print('Warning: ' + hl('no gate truthtables') + ' found!')
            return
        for n in range(len(self.gates)):
            self.bound_env[n] = self.structure.gate_truthtables[self.gates[n].node][ttix]
            self.gates[n].env = self.bound_env[n]

    def set_initial_value(self, values, ttix):
        val = np.zeros([len(self.node_idx), self.wire_capacity])
        if ttix != -1:
            self.set_truthtable(ttix)
        for input in self.structure.inputs:
            # Place the input at the promoter activity position, i.e. 0
            val[self.node_idx[input], 0] = values[reverse_input_encoding[input]]
        self.solver.set_initial_value(val)

    def set_solver(self, solver_class):
        self.solver = solver_class(self, None, False)

    def solve(self, tol, max_iter):
        if self.solver is not None:
            (values, err, iter) = self.solver.solve(tol, max_iter)
            # currently supports only one output
            # Place the output at the RNA concentration position, i.e. 0
            return (values[self.node_idx[list(self.structure.outputs)[0]], 1], err, iter)
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
                if n['type'].startswith('INPUT'):
                    self.stats['n_inputs'] += 1
                    self.inputs.add(n['id'])
                elif n['type'].startswith('OUTPUT'):
                    self.stats['n_outputs'] += 1
                    self.outputs.add(n['id'])
                else:
                    self.stats['n_nodes'] += 1
                    self.gates.add(n['id'])
            self.nodes = self.inputs | self.gates | self.outputs
            self.gate_truthtables = None
            if 'gate_truthtables' in circ:
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
            node_list = list(self.nodes)
            visited = dict(zip(list(self.nodes), [False for _ in range(len(self.nodes))]))
            while len(node_list) > 0:
                nix = 0
                for node in node_list:
                    good = True
                    for target in self.adjacency['out'][node]:
                        if not visited[target]:
                            good = False
                            break
                    if not good:
                        nix += 1
                        continue
                    self.order.insert(0, node)
                    visited[node] = True
                    node_list.pop(nix)
                    break
            if DEBUG_LEVEL > 1:
                print('Combinational order:')
                print(' -> '.join(list(map(hl, self.order))))
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
        self.dummys = []
        for k, v in assi.items():
            self.map_gtod[k] = v
            self.map_dtog[v] = k
        self.devices = set(self.map_dtog.keys())
        self.gates = set(self.map_dtog.values())

class library:
    # represents a connected gate lib
    def __init__(self, jstr):
        self.c = dict()
        libstr = json.loads(jstr)
        for entry in libstr:
            self.c[entry['class']] = entry['members']
    # return the JSON representation of this lib
    def __str__(self):
        content = list()
        for c, m in self.c.items():
            content += [dict()]
            content[-1]['class'] = c
            content[-1]['members'] = m
        return json.dumps(content)

# PRIVATE. Not intended to use manually
# unified base gate
class _base_gate:
    # a base class for a logic gate
    def __init__(self, node, dev):
        self.node = node
        self.dev = dev
        self.type = -2
        self.env = 0
        self.inputs = list()

# PRIVATE. Not intended to use manually
# The implicit OR gate is just a pseudo-gate, which outputs the YFP RPU
class _feedback_implicit_or_gate(_base_gate):
    # represents a NOR gate in the circuit
    def __init__(self, node, dev, parameters):
        self.node = node
        self.dev = dev
        self.type = 1
        self.mu = parameters['mu']
        self.gamma = parameters['gamma']
        self.delta = parameters['delta']
        #self.delta = 1.
        print(parameters)
    def out(self, input):
        # input contains three values
        # input[0] = sum_inwired activity
        # input[1] = sum_inwired RNA
        # input[2] = sum_outwired RNA * delta_outwired
        # output contains two values
        # output[0] = updated activity (nonexistant/unused)
        # output[1] = updated RNA
        # output[2] = updated RNA times delta
        output = np.zeros(3)
        output[1] = (self.mu*input[0])/(self.gamma + self.delta*input[1] + input[2])
        output[2] = self.delta*output[1]
        return output
    def __str__(self):
        return self.name + ': None (implicit OR)'

# The NOR gate with feedback is the gate proposed at IWBDA 2023
class _feedback_nor_gate(_base_gate):
    # represents a NOR gate in the circuit
    def __init__(self, node, dev, parameters):
        self.node = node
        self.dev = dev
        self.mu = parameters['mu']
        self.gamma = parameters['gamma']
        self.kappa = parameters['kappa']
        self.delta = parameters['delta']
        self.n = parameters['n']
        self.beta = parameters['beta']
        self.type = 0
        print(parameters)
    def out(self, input):
        # input contains three values
        # input[0] = sum_inwired activity
        # input[1] = sum_inwired RNA
        # input[2] = sum_outwired RNA * delta_outwired
        # output contains two values
        # output[0] = updated activity
        # output[1] = updated RNA
        # output[2] = updated RNA times delta
        output = np.zeros(3)
        output[1] = (self.mu*input[0])/(self.gamma + self.delta*input[1] + input[2])
        output[0] = (np.power(self.kappa, self.n)*(1. - self.beta))/(np.power(self.kappa, self.n) + np.power(output[1], self.n)) + self.beta
        output[2] = self.delta*output[1]
        return output
    def __str__(self):
        return self.name + ': ' + str(self.e)

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
    def solve(self, tol=10**(-2), max_iter=100):
        # Iterate the points using Banach's fixed point theorem
        # to estimate the error, compare old and new node values
        err = np.inf
        run = 0
        function = self.circuit.propagate
        vs = self.values
        while err > tol:
            #if True:#DEBUG_LEVEL > 1:
                #self._debug_step(vs, run, err)
            new_vs = function(vs)
            # error and max_iter updates
            err = np.max(np.abs(new_vs - vs))
            vs = new_vs
            run += 1
            if run >= max_iter:
                break
        if True:
            self._debug_step(vs, run, err)
        self.values = vs
        return (vs, err, run)
    def _debug_step(self, vs, run, err):
        print('---------\n' + head('Outputs') + ' (run = ' + hl(str(run)) + ', err < ' + hl(str(err)) + ', logic = ' + hl(str(self.circuit.bound_env[self.circuit.node_idx[list(self.circuit.structure.outputs)[0]]])) + '): ')
        for k, v in self.circuit.node_idx.items():
            if self.circuit.gates[self.circuit.node_idx[k]].type == 0:
                print('{k:<30}'.format(k = hl(k) + ': ' + hl(str(vs[v]))))
            elif self.circuit.gates[self.circuit.node_idx[k]].type == -2:
                print('{k:<30}'.format(k = hl(k) + ' (input)') + ': ' + hl(str(vs[v])))
            elif self.circuit.gates[self.circuit.node_idx[k]].type == 1:
                print('{k:<30}'.format(k = hl(k) + ' (implicit or)') + ': ' + hl(str(vs[v])))

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

class _directed_acyclic_graph:
    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges
    def is_valid(self):
        pass

def _list_intersect(a, b):
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

# JIT friendly version of integer power
def _power(bases, exponents):
    for n in range(len(bases)):
        b = bases[n]
        e = exponents[n] - 1
        for m in range(e):
            bases[n] *= b
