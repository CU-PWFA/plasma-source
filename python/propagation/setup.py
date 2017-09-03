#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 15:09:18 2017

@author: robert

Run from the propagation folder with the command:
    'python setup.py build_ext --inplace'
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

ext_modules=[
    Extension('laserv2',
              sources=['laserv2.pyx'],
              libraries=["m"],
              #extra_compile_args = ['-03', '-march=native', '-fopenmp'],
              #extra_link_args=['-fopenmp']
    )
]

setup(ext_modules=cythonize(ext_modules))