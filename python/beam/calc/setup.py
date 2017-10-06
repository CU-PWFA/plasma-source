#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 14:46:39 2017

@author: robert

Run from the python folder with the command:
    'python beam/calc/setup.py build_ext --inplace'
"""

import numpy
import venv
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


ext_modules=[
    Extension('beam.calc.*',
              sources=['beam/calc/*.pyx'], #Compile entire module
              #sources=['beam/calc/ionization.pyx'],
              #sources=['beam/calc/laser.pyx'], #Compile specific files
              libraries=["m", "fftw3"],
              include_dirs=[numpy.get_include(), 
                            os.path.join(venv.sys.base_prefix, 'include')],
              extra_compile_args = ['-march=native', '-fopenmp', '-O3',
                                    '-lfftw3'],
              extra_link_args=['-fopenmp'],
              define_macros=[('CYTHON_TRACE', '1')]
    )
]

setup(ext_modules=cythonize(ext_modules))
