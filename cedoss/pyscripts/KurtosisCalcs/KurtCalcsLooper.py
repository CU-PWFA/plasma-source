#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 12:05:37 2018

Calculate the kurtosis for given functions
This one has more loops and less plots

@author: chris
"""

import numpy as np
import KurtCalcs
import matplotlib.pyplot as plt

Jmax = 1.5e-8
em = 5e-6/2e4
#Bm = 1.1
num = int(1e5)

Bm_arr = np.linspace(1.0,2.0,20)
ave_kurt = np.zeros(len(Bm_arr))
for j in range(len(Bm_arr)):
    Bm = Bm_arr[j]
    kurt_arr = np.zeros(10)
    for i in range(len(kurt_arr)):
        kurt_arr[i] = KurtCalcs.CalcKurt(KurtCalcs.FullDist, Jmax, em, Bm, sample_num = num)
        print(kurt_arr[i])
    print(Bm)
    print(" Ave: ",np.average(kurt_arr))
    print()
    ave_kurt[j] = np.average(kurt_arr)
    

plt.plot(Bm_arr, ave_kurt)
plt.title("Average (N=10) kurtosis vs. B-mag")
plt.xlabel("B-mag")
plt.ylabel(r'$\kappa$'+" from Raubenheimer, 1996, Eqn. 21")
plt.grid(); plt.show()