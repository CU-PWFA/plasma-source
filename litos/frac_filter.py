#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 17:28:57 2017

@author: litos
"""

import numpy as np

def frac_filter(x,frac=1.0,medmean=''):
        
    if medmean.lower() == 'mean':
        dx = x - np.mean(x)
    elif medmean.lower() == 'med':
        dx = x - np.median(x)
    else :
        dx = x - (np.mean(x)+np.median(x))/2
    
    ifilt = sorted(range(len(dx)),key=lambda ifilt: dx[ifilt])
    nx = round(frac*len(ifilt))
    ifilt = ifilt[:nx]
     
    return ifilt