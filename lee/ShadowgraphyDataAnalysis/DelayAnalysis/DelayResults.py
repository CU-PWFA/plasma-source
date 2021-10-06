#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 17:41:00 2021

@author: valentinalee
"""


import numpy as np
import matplotlib.pyplot as plt

#%%
#[Delay, PeakToTrough(P-V)/V, PeakWidth]
Results= np.array([[90.5, 1.39526, 45.812], [88, 1.3045, 39.8], [80.5, 0.6603972, 40.144], \
          [73, 3.437287, 43.428], [63, 2.4684396, 38.656], [53, 6.0131557, 40.504], \
          [43, 7.246088, 43.248], [33, 8.163498, 40.0], [23, 2.619355, 37.96], \
          [15.5, 4.707629, 41.284], [13, 2.6043, 39.928], [-6, 2.847013, 41.604], \
          [-18, 2.6738628, 41.144], [-42, 2.923777, 43.488], [-54, 4.2348, 41.82], \
          [30.5, 1.83114, 40.108]])

#[-30, 8.09482, 44.012], 
#%%
plt.figure(1)
#plt.plot((Results[:, 0])*-1, Results[:, 1], 'o')
plt.plot((Results[:, 0]-90.5)*-1*4/30, Results[:, 1], 'o')
plt.xlabel('Delay (ns)')
plt.ylabel('(Peak-Trough)/Trough')
#%%
plt.figure(2)
plt.plot((Results[:, 0]-90.5)*-1*4/30, Results[:, 2], 'o')