from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Four State Promoter Model',
    ext_modules=cythonize("four_state_promoter_model.py", force=True),
)


# To run setup
# python setup.py build_ext --inplace