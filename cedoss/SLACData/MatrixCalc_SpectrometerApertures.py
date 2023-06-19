#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 17:33:07 2023

Trying to plot out the beam envelope as it propagates through the spectrometer.

Just an early version for general spot sizes.  A later implementation into
WARGSim should also try to calculate the beam loss.

Sorry, finished the thesis and then never got around to this :(
I'll at least leave the spectrometer apertures from Doug S. somewhere

@author: chris
"""

import numpy as np
import matplotlib.pyplot as ply

import MatrixCalc

###Define all the constants

###Define the zones of the spectrometer, and fully plan out the total steps

###In each zone, propagate a beam using the step size

###Plot out the sigma, and plot out the apertures.