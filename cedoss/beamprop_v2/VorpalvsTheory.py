#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 11:53:41 2018

Plots sim data against theory

@author: chris
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import BeamPropFuncs as PProp
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import CalcEmitGrowth as W2

e0_arr = np.array([5.301,  5.301,  5.301,  5.301])
ef_arr = np.array([6.629, 17.629, 33.882, 50.481])
thick_arr = np.array([50, 200, 400, 600])

tpl_n = 10. #e17 cm-3
betastar = 0.05
gammab = PProp.def_gamma

bmag_arr = ef_arr/e0_arr
len_arr = np.linspace(20, 700, 100)
k = Foc.Calc_K(tpl_n*1e17, gammab)*100*100
sigmaE = 0.005

bmag_w2_arr_thick = np.zeros(len(len_arr))
for x in range(len(len_arr)):
    w2 = W2.ThickW2_UnNormalized(k, len_arr[x]*1e-6, betastar, 0)
    bmag_w2_arr_thick[x] = W2.CalcEmit(w2, sigmaE)
    
plt.plot(thick_arr, bmag_arr)
plt.plot(len_arr, bmag_w2_arr_thick)