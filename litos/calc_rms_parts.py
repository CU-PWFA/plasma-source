#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:11:35 2017

@author: litos
"""

import numpy as np
from calc_rms import calc_rms

def calc_rms_parts(parts,frac=1.0):
    
    npart = len(parts)
    x  = np.zeros(npart)
    xp = np.zeros(npart)
    y  = np.zeros(npart)
    yp = np.zeros(npart)
    z  = np.zeros(npart)
    gb = np.zeros(npart)
    for i in range(0,npart):
        [x[i],xp[i],y[i],yp[i],z[i],gb[i]] = parts[i][:]
    
    rms_x  = calc_rms(x,frac)
    rms_xp = calc_rms(xp,frac)
    rms_y  = calc_rms(y,frac)
    rms_yp = calc_rms(yp,frac)
    rms_z  = calc_rms(z,frac)
    rms_gb = calc_rms(gb,frac)
     
    return [rms_x,rms_xp,rms_y,rms_yp,rms_z,rms_gb]