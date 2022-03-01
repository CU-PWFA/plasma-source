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
superpath = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/'

"""
#path = '/home/chris/Desktop/10e16/'
path = '/home/chris/Desktop/NERSC_Deflection_Test/'
path = '/home/chris/Desktop/NERSC_Deflection_July/'
npcase = 2e16
multfac = 2
vmax = 1e18*(npcase/3.7e17)*multfac
vmin = 3e16*(npcase/3.7e17)/multfac/7
params = {'drive' : 'rhoDrive',
          'witness' : 'rhoWitness',
          'plasma' : 'electrons',
          'dumpInd' : 10,
          'path' : path,
          #'simName' : 'ThinPlasmaLens3D',
          'simName' : 'PTPL_Gradient',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'vmax' : vmax,
          'vmin' : vmin
          }
plot.drive_witness_density_new(params)
"""

#path = '/home/chris/Desktop/WakeShape_LinGrad/'
#path = superpath + 'NERSC_Sep_Grad/'
#path = '/home/chris/Desktop/NERSC_LongiFieldFix1/'
#path = '/home/chris/Desktop/NERSC_Deflection_July/'
path =  '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_n2e16_g8e17/'
#path =  '/media/chris/New Volume/VSimRuns/AugustLinearGradient/Tests/NERSC_n3e17_g0/'
#path =  '/media/chris/New Volume/VSimRuns/AugustLinearGradient/NERSC_n1e17_g0/'
#path = superpath + 'NERSC_Dec_Grad/'
params = {'drive' : 'rhoDrive',
          'witness' : 'rhoWitness',
          'plasma' : 'electrons',
          'path' : path,
          'dumpInd' : 5,
          #'simName' : 'ThinPlasmaLens3D',
          'simName' : 'MatchedBeams',
          #'simName' : 'PTPL_Gradient',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05
          }
#plot.drive_witness_density_new(params)
plot.drive_witness_density_paperplot(params)

"""
path = '/home/chris/Desktop/SimulationRepo/TestCases/3D/SingleBunch/5e17/200um/'
path = '/home/chris/Desktop/NERSC_LIN_Aug/'
path = '/media/chris/New Volume/VSimRuns/NERSCResults/Paper1/FACETII_NERSC_Run4/'

params = {'drive' : 'rhoDrive',
          'plasma' : 'electrons',
          'dumpInd' : 3,
          'path' : path,
          #'simName' : 'MatchedBeams',
          'simName' : 'PTPLDoubleTanh',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05
          }
plot.drive_witness_density_new(params)
"""
"""
path = '/media/chris/New Volume/VSimRuns/NERSCResults/3D/TwoBunch/1e18/50um/'
path = '/home/chris/Desktop/thick4/'
params = {'drive' : 'rhoDrive',
          'witness' : 'rhoWitness',
          'plasma' : 'electrons',
          'dumpInd' : 10,
          'path' : path,
          #'simName' : 'ThinPlasmaLens3D',
          'simName' : 'PTPLDoubleTanh',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05
          }
plot.drive_witness_density_new(params)
"""
"""
path = '/home/chris/Desktop/SFQED_NERSC_VaryDen/5e15/'
#path = '/home/chris/Desktop/SimulationRepo/emittance_preservation/simulations/PTPL_FacetIIBeams/'
params = {'drive' : 'rhoDrive',
          'witness' : 'rhoWitness',
          'plasma' : 'electrons',
          'dumpInd' : 5,
          'path' : path,
          'simName' : 'PTPLDoubleTanh',
          'zoom' : 4.0,
          'alphaCutoff' : 0.001
          }
plot.drive_witness_density_new(params)

"""