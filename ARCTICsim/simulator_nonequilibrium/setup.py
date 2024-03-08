from setuptools import setup
from Cython.Build import cythonize

SIM_DIR = "ARCTICsim/simulator_noneqilibrium/"
SIM_DIR = ""

setup(
    name='Circuit Evaluator',
    ext_modules=cythonize([SIM_DIR + file for file in ["simulator/circuit_evaluator_new.py",
                                                       "simulator/circuit.py",
                                                       "models/*.py"]], force=True),
)

# To run setup
# python setup.py build_ext --inplace


# cythonizing the circuit evaluator does not provide any speed up.
