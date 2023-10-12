from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Circuit Evaluator',
    ext_modules=cythonize("circuit_evaluator.py", force=True),
)


# To run setup
# python setup.py build_ext --inplace


# cythonizing the circuit evaluator does not provide any speed up.