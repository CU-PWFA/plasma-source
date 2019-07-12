#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 13:27:34 2019

This version assumes we know the distance between the lens and the new waist,
more useful if you know there is a minimum distance between the lens and 
the plasma

Specifically, this had to assume an initial beta waist.  Another version
could assume that d=0 to solve for what bi could be.  The mathematica code
to go along with this is RestrictedTPL_solve_d_L.nb

@author: chris
"""

import numpy as np
import sys
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

k = 9047.683978841857
bf = 0.04214
bi = 0.10
zw = 0.7041

print("First solution")
l = (bi*k*zw-np.sqrt(bf*bi*k**2*(bf**2-bf*bi+zw**2)))/(bi*k**2*(bf**2+zw**2))
print("L",l*1e6,"um")
d = -(np.sqrt(bf*bi*k**2*(bf**2-bf*bi+zw**2))/(bf*k))
print("d",d,"m")

print("")
print("Second solution")
l = (bi*k*zw+np.sqrt(bf*bi*k**2*(bf**2-bf*bi+zw**2)))/(bi*k**2*(bf**2+zw**2))
print("L",l*1e6,"um")
d = (np.sqrt(bf*bi*k**2*(bf**2-bf*bi+zw**2))/(bf*k))
print("d",d,"m")