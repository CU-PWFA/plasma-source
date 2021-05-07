#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 12:31:39 2019

@author: valentina_lee
"""

#%%
import numpy as np;

#%%
Pinmbar=14.2
PinPa= Pinmbar*100
kB=1.3807e-23
T=293
NoverV= PinPa/kB/T
DenIncm3=NoverV/(100**3)

print(DenIncm3)