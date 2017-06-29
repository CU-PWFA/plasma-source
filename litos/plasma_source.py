#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 15:18:15 2017

@author: mike
"""

import numpy as np
from collections import defaultdict

def calc_dgds(dgds0,npl0,npl,model=0):
    if model==0:
        # wake strength ~sqrt(np), phase ~sqrt(np)
        dgds = dgds0*np.sqrt(npl/npl0)*(2*np.sqrt(npl/npl0)-1)
    else:
        dgds = dgds0
        
    return dgds

def make_ramp(s,updn,shape,hw,s_top,npl0,dgds0):
    ramp = defaultdict(dict)
    
    ramp["L"]     = np.abs(s[-1]-s[0])
    ramp["updn"]  = updn
    ramp["shape"] = shape
    ramp["hw"]    = hw
    ramp["s_top"] = s_top
    ramp["npl0"]  = npl0
    ramp["dgds0"] = dgds0
    
    # normal gauss. func.
    if shape.lower() == 'gauss':
        sig = hw/(np.sqrt(2*np.log(2)))
        # up-ramp
        if updn.lower() == 'up':
            npl = npl0*((s<s_top)*np.exp(-((s-s_top)**2)/(2*(sig**2))) +\
                       (s>=s_top)*1)
        # down-ramp
        else:
            npl = npl0*((s>s_top)*np.exp(-((s-s_top)**2)/(2*(sig**2))) +\
                       (s<=s_top)*1)
    # sigmoidal func.
    elif shape.lower() == 'sigmoid':
        sig = hw
        # up-ramp
        if updn.lower() == 'up':
            npl = npl0*(1/(1+np.exp(-(s-(s_top-2*sig))/(sig/4))))
        # down-ramp
        else:
            npl = npl0*(1/(1+np.exp(+(s-(s_top+2*sig))/(sig/4))))
    # no ramp
    else:
        npl = np.zeros(len(s))

    ramp["s"]    = s
    ramp["npl"]  = npl
    ramp["dgds"] = calc_dgds(dgds0,npl0,npl)

    return ramp

def make_bulk(s,npl0,dgds0):
    bulk = defaultdict(dict)
    
    bulk["L"]     = np.abs(s[-1]-s[0])
    bulk["npl0"]  = npl0
    bulk["dgds0"] = dgds0
    bulk["s"]     = s
    bulk["npl"]   = npl0*np.ones(len(s))
    bulk["dgds"]  = dgds0*np.ones(len(s))
    
    return bulk

def make_plasma(bulk,up_ramp=0,dn_ramp=0):
    plasma = defaultdict(dict)
    plasma["s"]    = []
    plasma["npl"]  = []
    plasma["dgds"] = []
    
    if up_ramp!=0:
        plasma["up_ramp"] = up_ramp
        plasma["s"]       = up_ramp["s"]
        plasma["npl"]     = up_ramp["npl"]
        plasma["dgds"]    = up_ramp["dgds"]
    
    plasma["bulk"] = bulk
    plasma["s"]    = np.append(plasma["s"][:-1],\
                                  plasma["s"][-1]+bulk["s"])
    plasma["npl"]  = np.append(plasma["npl"][:-1],bulk["npl"])
    plasma["dgds"] = np.append(plasma["dgds"][:-1],bulk["dgds"])
    
    if dn_ramp!=0:
        plasma["dn_ramp"] = dn_ramp
        plasma["s"]       = np.append(plasma["s"][:-1],\
                                  plasma["s"][-1]+dn_ramp["s"])
        plasma["npl"]     = np.append(plasma["npl"][:-1],dn_ramp["npl"])
        plasma["dgds"]    = np.append(plasma["dgds"][:-1],dn_ramp["dgds"])

    return plasma

