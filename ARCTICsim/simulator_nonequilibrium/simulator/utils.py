"""
Author: Erik Kubaczka
"""

import json
import sys
from configparser import ConfigParser

from ARCTICsim.simulator_nonequilibrium.simulator.libsim_coop import json_file


#####################################################################################################
#                                                                                                   #
# The subsequent classes represent extensions to classes in libsim_coop.py from Nicolai Engelmann   #
#                                                                                                   #
#####################################################################################################
class JsonFile(json_file):
    def __init__(self, path=None, content=None, load=True):
        super().__init__(path, read=load if content is None else False)

        if path is None:
            self.content = content

        self.data = None

        if load:
            self.load()

    def load(self):
        self.data = json.loads(self.content)

    def dump(self, data):
        if self.path is None or self.path == "":
            raise Exception("No path specified to dump to json file.")
        self.content = json.dumps(data)
        self.data = data
        self.write()


#####################################################################################################
#                                                                                                   #
# The subsequent is taken from Nicolai Engelmann's thermodynamic simulator                          #
#                                                                                                   #
#####################################################################################################

def load_settings(path):
    s = ConfigParser(allow_no_value=True)
    s.read(path)
    s = dict({k: v for sec in s.sections() for k, v in s[sec].items()})
    c = dict()
    t = dict()
    for k, v in s.items():
        s[k], t[k], c[k] = tuple(map(str.strip, v.split('#', 2)))
        s[k] = type_dict[t[k]](s[k])
    return s, t, c


type_dict = dict({
    'int': int,
    'str': str,
    'bool': bool,
    'float': float,
})


class communication_wrapper:
    def __init__(self, input_file, output_file, prefix=None):
        self.i = None
        self.o = None
        self.prefix = prefix
        if input_file is None:
            self.i = sys.stdin
        else:
            pass  # stored and temporary files not supported yet
        if output_file is None:
            self.o = sys.stdout
        else:
            pass  # stored and temporary files not supported yet

    def readline(self):
        # first signal, that you are ready, by printing the prefix
        if self.prefix is not None:
            self.o.write(self.prefix + ' ')
            self.o.flush()
        return self.i.readline()

    def writeline(self, line):
        self.o.write(line + '\n')
