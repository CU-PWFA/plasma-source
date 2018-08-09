#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 12:10:26 2018

Calculates things related to W2 for plasma lenses

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

focal = 0.001
kl = 1/focal

bstar = 0.10
b = bstar * kl

dsep = -0.10
d = dsep * kl

sigmaE = 0.0058

#W2 = (b**2+d**2)**2/b**2
W2 = 19666.5

bmag = 1 + 1/2*W2*sigmaE**2

print(bmag)