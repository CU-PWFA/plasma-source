#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:59:08 2018

To get all the good stuff, get 'electron' 'rhoWitness' 'rhoDrive for all the dumps
then just the final 'WitnessBeam' and 'DriveBeam' dump.

@author: chris
"""

import sys
sys.path.insert(0, "../../python")
import os
import numpy as np
from vsim import plot
import matplotlib.pyplot as plt
"""
path = '/home/chris/Desktop/SimulationRepo/TestCases/3D/SingleBunch/5e17/200um/'
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
"""
path = '/media/chris/New Volume/VSimRuns/NERSCResults/3D/TwoBunch/1e18/50um/'
params = {'drive' : 'rhoDrive',
          'witness' : 'rhoWitness',
          'plasma' : 'electrons',
          'dumpInd' : 3,
          'path' : path,
          'simName' : 'ThinPlasmaLens3D',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05
          }
plot.drive_witness_density_new(params)
"""
#"""
path = '/home/chris/Desktop/FACETII_NERSC_Run3/'
#path = '/home/chris/Desktop/SimulationRepo/emittance_preservation/simulations/PTPL_FacetIIBeams/'
params = {'drive' : 'rhoDrive',
          'witness' : 'rhoWitness',
          'plasma' : 'electrons',
          'dumpInd' : 8,
          'path' : path,
          'simName' : 'PTPLDoubleTanh',
          'zoom' : 4.0,
          'alphaCutoff' : 0.01
          }
plot.drive_witness_density_new(params)

#"""