#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:48:04 2017

Functions for creating a 1D/2D elliptical tanh to represent a thin plasma lens.

For 1D, pass the distance values through DoubleTanh along with parameters for
the tanh profile.

For 2D, first pass the the 2D distances through EllipDist along with the scaling
parameter for the ellipse, then pass the returned 'effective distance' into
DoubleTanh just as you would for 1D.

DefaultParams contains tanh parameters from the fit of the current plasma lens
I am working with

An example of this distribution can be seen by running the script
"cedoss/pyscripts/ApproxDensity.py"

@author: chris
"""

import numpy as np

#Generic function for a double tanh profile.  Flat in the center with ramps
# on either side which go to zero.  Centered around x = 0
# The parameters 'a' and 'b' can be roughly thought of as the flat top and ramp
# length, repsectively; but in reality 'a' is the Half Width at Half Max and
# 'b' controls the slope at x = +/- a.  As b ---> 0 the tanh profile becomes
# a square profile.  'b' is usually an order of magnitude smaller than 'a'
#
#  a - Half Width Half Max
#  b - ~ramp length
#  n_0 - density at center
#  x - array of distances
def DoubleTanh(x, a, b, n_0):
    return ((.5 + .5*np.tanh((x + a)/b)) * 
            (.5 - .5*np.tanh((x - a)/b))) * n_0
            
#Returns the effective distance for an ellipse, y^2 * (scl*z)^2
#
#  x,y - array of distances in the elliptical plane.
#  scl - scaling factor between the narrow and wide waists
#        typically if x is the narrow waist and y is the wide waist
#        then scl is roughly the ratio of waist_x / waist_y
def EllipDist(x, y, scl):
    return np.sqrt(np.square(x) + np.square(scl * y))
   
#Some example values for elliptical tanh's I find with the current setup. These
# particular values correspond to a laser pulse propagating through uniform
# density of 1e17 cm^-3.  The plasma can be assumed uniform in the direction
# orthogonal to the elliptical plane, as in this direction the plasma is
# about 100x wider than the narrow waist due to the Rayleigh length
#
#  See functions above for descriptions and usage of these parameters
def DefaultParams():
    a = 13.3     #In units microns
    b = 2.14     #In units microns
    n_0 = 1e17   #In units cm^-3
    scl = 0.0587 #Dimensionless
    return [a, b, n_0, scl]