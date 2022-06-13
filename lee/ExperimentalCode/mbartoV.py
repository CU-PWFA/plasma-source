#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 10:34:38 2021

@author: valentinalee
"""

import numpy as np
#%%
mv_arr = np.array([30,42,55,66])
mbar_arr = np.array([1,2,3,4])
k = np.polyfit(mbar_arr, mv_arr ,1)
print(k)
mbar_exp = 8
mv_exp = k[0]*mbar_exp + k[1]
print(mv_exp*1e-3, 'V')
