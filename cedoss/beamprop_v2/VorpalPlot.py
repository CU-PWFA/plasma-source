#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:59:08 2018

@author: chris
"""

import sys
sys.path.insert(0, "../../python")
import os
import numpy as np
from vsim import plot
import matplotlib.pyplot as plt
"""
path = '/home/chris/Desktop/SimulationRepo/TestCases/3D/200um/'
params = {'drive' : 'rhoDrive',
          'plasma' : 'electrons',
          'dumpInd' : 4,
          'path' : path,
          'simName' : 'ThinPlasmaLens3D_Single',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05
          }
plot.single_drive_density(params)
"""

path = '/home/chris/Desktop/SimulationRepo/TestCases/3D/TwoBunch/1e18/200um/'
params = {'drive' : 'rhoDrive',
          'witness' : 'rhoWitness',
          'plasma' : 'electrons',
          'dumpInd' : 4,
          'path' : path,
          'simName' : 'ThinPlasmaLens3D',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05
          }
plot.drive_witness_density(params)