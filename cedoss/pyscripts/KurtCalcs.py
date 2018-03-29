#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 12:05:37 2018

Calculate the kurtosis for given functions

@author: chris
"""

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import random
import scipy.special as spec

case = 1

def I_0(arg):
    return spec.i0(arg)
    """
    for i in range(max_order+1):
        res = res + (1/4*arg**2)**i/(math.factorial(i))**2
    """

def BaseDist(J, em, Bm):
    return 1/em*np.exp(-Bm*J/em)

def BesselTerm(J, em, Bm):
    return I_0(J/em*np.sqrt(Bm**2-1))

def FullDist(J, em, Bm):
    return BaseDist(J,em,Bm)*BesselTerm(J,em,Bm)

def CalcKurt(func, Jmax, em, Bm, sample_num = int(1e5), plot=False):
    Jmin = 0
    
    set_arr = np.linspace(0,1,num=int(1e5))
    Jset_arr = np.linspace(0,1,num=len(set_arr))
    
    ymin = func(Jmax, em, Bm)
    ymax = func(Jmin, em, Bm)
    """
    numSteps = 1000
    for i in range(numSteps):
        J = Jmin + (Jmax - Jmin) * float(i) / numSteps
        y = BaseDist(J, em, Bm)
        if y < ymin: ymin = y; print("hi")
        if y > ymax: ymax = y; print("hello")
    """
    
    for i in range(len(set_arr)):
        if case == 0:
            x = random.gauss(0, 5e-6)
            xp = random.gauss(0, 5e-6)
            Jset_arr[i] = (x**2+xp**2)/2
        if case == 1:
            while True:
                # generate a random number between 0 to 1
                Jr = random.random()
                yr = random.random()
                J = Jmin + (Jmax - Jmin) * Jr
                y = ymin + (ymax - ymin) * yr
                if y <= func(J, em, Bm):
                    Jset_arr[i] = J
                    break
    
    kurt = stats.kurtosis(Jset_arr,0,False,True)
    if plot:
        plt.hist(Jset_arr,20)
        J_arr = np.linspace(0,Jmax,1000)
        base_arr = func(J_arr, em, Bm)
        plt.plot(J_arr,base_arr/base_arr[0]*sample_num/2)
        plt.show()
    return kurt
