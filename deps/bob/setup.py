from setuptools import setup
from Cython.Build import cythonize

setup(
    name='bob',
    ext_modules=cythonize("src/bob.pyx")
)