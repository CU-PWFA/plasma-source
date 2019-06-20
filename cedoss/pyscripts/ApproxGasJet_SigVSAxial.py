#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 14:00:09 2019

I felt bad for not including sig varying with axial distance.  To rectify this
I simply measured three values of sigma in the domains of 5mm +/- 1.2mm and
will do a quad fit.

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

"""This was for the WedgeRAS R2 simulation
y_arr = [-1200, 0, 1200] #um
sig_arr = [2858.44472812, 3536.31577167, 4187.60231364] #um
"""
y_arr = [-1200, 0, 1200] #um
sig_arr = [2319.33640609, 2986.91262632, 3639.82876814] #um

pfit = np.polyfit(y_arr, sig_arr, 2)
print('sig(y)=',pfit[2],'+',pfit[1],'y','+',pfit[0],'y**2')

plt.plot(y_arr, sig_arr)
plt.show()