#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:20:26 2018

Runs KurtCalc once and plots.  Does other misc things with KurtCalcs

@author: chris
"""

import KurtCalcs
import matplotlib.pyplot as plt
import numpy as np

Jmax = 5.0e-9
em = 5e-6/2e4
Bm = 2.0
num = int(1e5)
"""
kurt = KurtCalcs.CalcKurt(KurtCalcs.FullDist, Jmax, em, Bm, sample_num = num, plot=True)
print(kurt)
"""
jmax_arr = np.linspace(5e-9,5e-8, num=10)
kurt_arr = np.zeros(len(jmax_arr))
for i in range(len(jmax_arr)):
    kurt_arr[i] = KurtCalcs.CalcKurt(KurtCalcs.FullDist, jmax_arr[i], em, Bm, sample_num = num, plot=False)
    print(kurt_arr[i])

plt.plot(jmax_arr,kurt_arr)
plt.title("Convergence: Kurtosis vs Jmax w/ Bm = "+str(Bm))
plt.xlabel("Jmax")
plt.ylabel("Kurtosis of one distribution")
plt.grid(); plt.show()


"""
J_arr = np.linspace(0,Jmax,100)
bess_arr = KurtCalcs.BesselTerm(J_arr, em, Bm)
plt.plot(J_arr, bess_arr)
plt.title("I_0 for B-mag = " + str(Bm))
plt.xlabel("J")
plt.ylabel("I_0")
plt.grid(); plt.show()
"""