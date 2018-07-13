#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:59:52 2018

My implementation of Robert's code to plot phase space of VSim species

@author: chris
"""

import sys
import os
sys.path.insert(0,"../../python")
from vsim import plot

if __name__ == "__main__":
    #arg = 'Drive_Witness_Ramps_WitnessBeam_2.h5'
    arg = sys.argv[1]
    
    strlist = arg.split("_")
    dumpInd = int(strlist[-1].split(".")[0])
    species = strlist[-2]
    simname = "_".join(strlist[:-2])
    
    #cwd = '/home/chris/Desktop/'
    cwd = os.getcwd()
    
    params = {
        'species' : species,
        'dumpInd' : dumpInd,
        'path' : cwd,
        'simName' : simname,
        'cutoff' : 0.0
        }
    print(params)
    #plot.phase_space(params)