#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 11:20:11 2020

quick plot of how ellip scaling is for eccentricity.

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

rx = 1.0
ry_arr = np.linspace(1.0,2.0,250)
#delt_arr = np.zeros(len(ry_arr))
#rmux_arr = np.zeros(len(ry_arr))

delt_arr = 0.5 + 0.5*ry_arr/rx
rmux_arr = 2*rx*ry_arr/(rx+ry_arr)

plt.title("Factor comp.:  Mike's Analytical (blue) VS Chris's Empirical (red)")
plt.plot(ry_arr/rx, delt_arr, c='b', label = r'$1+\delta/2$')
plt.plot(ry_arr/rx, rmux_arr, c='r', label = r'$r_\mu / r_x$')
plt.ylabel(r'$f(r_y/r_x)$')
plt.xlabel(r'$r_y/r_x$')
plt.grid(); plt.legend(); plt.show()