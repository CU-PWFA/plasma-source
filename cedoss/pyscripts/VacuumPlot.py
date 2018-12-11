#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 13:33:03 2018

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

minarr = np.array([0,2,5,10,15,17,20,25,30,35,40,45,55,70])
prsarr = np.array([3.15e3, 3.75e1, 6.75e-1, 1.06e-1, 1.03e-1, #up to 15 min
                   6.15e-4, 1.06e-4, 2.50e-5, 1.45e-5, 8.95e-6, #up to 35 min
                   6.71e-6, 5.32e-6, 4.50e-6, 4.10e-6])

plt.semilogy(minarr, prsarr)
plt.title("Vacuum Chamber Pressure vs Time")
plt.xlabel("Time [min]")
plt.ylabel("Pressure [mbar]")
plt.grid(); plt.show()