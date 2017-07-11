#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 13:44:27 2017

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

den = 5.73
top = 20
ramp = 10

#My initial guesses
"""
a = top/2 + ramp/2 #15.0
b = ramp/4         #2.5
"""

#From ThrDim.FitDataDoubleTanh
a = 211.81085591
b = 49.9021692099
den = 8.05244671013

window = 350

def f(x):
    first = .5 + .5*np.tanh((x+a)/b)
    second= .5 - .5*np.tanh((x-a)/b)
    return (first * second) * den

x = np.arange(-window,window,.01)
plt.plot(x,f(x))

"""
a = 15.1517144584
b = 2.59324260261
den = 8.05244671013

plt.plot(x,f(x))
"""