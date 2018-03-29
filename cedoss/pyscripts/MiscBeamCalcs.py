#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 13:23:12 2018

Calculates the effective size of a beam (sigma_r) away from waist

@author: chris
"""

import numpy as np

zv = -0.36
z = -0.40

sigv = 5.9e-6
en = 7.0e-6
gb = 19569.5

sig = sigv*np.sqrt(1+np.square(en/gb*(z-zv)/(sigv)**2))
print(sig)