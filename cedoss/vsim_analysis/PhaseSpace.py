#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:59:52 2018

My implementation of Robert's code to plot phase space of VSim species

@author: chris
"""

import sys
sys.path.insert(0,"../../python")
from vsim import plot

params = {
        'species' : 'WitnessBeam',
        'dumpInd' : 1,
        'path' : '/home/chris/Desktop/VSim_NERSC_Template_2dTest/',
        'simName' : 'Drive_Witness_Ramps',
        'cutoff' : 0.0
        }
plot.phase_space(params)