#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 15:09:18 2017

@author: robert

Run from the python folder with the command:
    'python propagation/setup.py build_ext --inplace'
"""

import numpy
import venv
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


ext_modules=[
    Extension('propagation.laserv2',
              sources=['propagation/laserv2.pyx'],
              libraries=["m", "fftw3"],
              include_dirs=[numpy.get_include(), 
                            os.path.join(venv.sys.base_prefix, 'include')],
              extra_compile_args = ['-march=native', '-fopenmp', '-O3',
                                    #'-DCYTHON_TRACE_NOGIL=1',
                                    #'-I/home/robert/anaconda3/envs/CU-PWFA/include', 
                                    '-lfftw3'],
              extra_link_args=['-fopenmp'],
              define_macros=[('CYTHON_TRACE', '1')]
    )
]

setup(ext_modules=cythonize(ext_modules))