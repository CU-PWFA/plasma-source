#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:18:05 2017

@author: litos
"""

import numpy as np
from draw_ellipse import draw_ellipse


eps   = 1e-6
beta  = 1
alpha = 0
gamma = (1+(alpha**2))/beta
T     = np.array([beta,alpha,gamma])
sig   = np.sqrt(eps*beta)

[x,xp] = draw_ellipse(eps,T)

print(x)
print(xp)