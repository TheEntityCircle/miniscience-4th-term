from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [

    Extension("using_cython", ["using_cython.pyx"]),

    Extension("using_cython_and_openmp", ["using_cython_and_openmp.pyx"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'])
]

setup(
    ext_modules = cythonize(extensions)
)
