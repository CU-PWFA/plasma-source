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

def make_ramp(s,updn,shape,hw,top_loc,npl0,dgds0):
    ramp = defaultdict(dict)
    
    ramp["L"]       = np.abs(s[-1]-s[0])
    ramp["updn"]    = updn
    ramp["shape"]   = shape
    ramp["hw"]      = hw
    ramp["top_loc"] = top_loc
    ramp["npl0"]    = npl0
    ramp["dgds0"]   = dgds0
    
    # check for zero-length ramp
    if hw==0:
        
        # up-ramp
        if updn.lower() == 'up':
            npl = npl0*((s<top_loc)*0 + (s>=top_loc))
        # down-ramp
        else:
            npl = npl0*((s>top_loc)*0 + (s<=top_loc))
    
    # non-zero-lenght ramps
    else:
        # normal gauss. func.
        if shape.lower() == 'gauss':
            # define gauss. sigma
            sig = hw/(np.sqrt(2*np.log(2)))            
            # up-ramp
            if updn.lower() == 'up':
                npl = npl0*((s<top_loc)*np.exp(-((s-top_loc)**2)/(2*(sig**2))) +\
                           (s>=top_loc)*1)
            # down-ramp
            else:
                npl = npl0*((s>top_loc)*np.exp(-((s-top_loc)**2)/(2*(sig**2))) +\
                           (s<=top_loc)*1)

        # sigmoidal func.
        elif shape.lower() == 'sigmoid':
            # up-ramp
            if updn.lower() == 'up':
                npl = npl0*(1/(1+np.exp(-(s-(top_loc-2*hw))/(hw/4))))
            # down-ramp
            else:
                npl = npl0*(1/(1+np.exp(+(s-(top_loc+2*hw))/(hw/4))))
                
        # trapezoidal func.
        elif shape.lower() == 'trap': # trapezoidal func.
            # up-ramp
            if updn.lower() == 'up':
                npl = npl0*((s<=top_loc-2*hw)*0 +\
                       (s>top_loc-2*hw)*(s<top_loc)*(s-(top_loc-2*hw))/(2*hw) +\
                       (s>=top_loc)*1)
            # down-ramp
            else:
                npl = npl0*((s>=top_loc+2*hw)*0 +\
                       (s<top_loc+2*hw)*(s>top_loc)*((top_loc+2*hw)-s)/(2*hw) +\
                       (s<=top_loc)*1)

        # gompertz function
        elif shape.lower() == 'gomp':
            # up-ramp
            if updn.lower() == 'up':
                npl = npl0*(1-np.exp(-np.exp((s-top_loc)/hw)))
            # down-ramp
            else:
                npl = npl0*(1-np.exp(-np.exp((top_loc-s)/hw)))

        # lorentzian function
        elif shape.lower() == 'lorentz':
            # up-ramp
            if updn.lower() == 'up':
                npl = npl0*(hw/np.pi)/((s-top_loc)**2+hw**2)
            # down-ramp
            else:
                npl = npl0*(hw/np.pi)/((top_loc-s)**2+hw**2)
                
        # Xu PRL 2016 ramp shape 3
        elif shape.lower() == 'xu3':
            # define sigma parameter
            sig = 2*hw
            # up-ramp
            if updn.lower() == 'up':
                npl = npl0*((s<top_loc)*(1/(1-2*(s-top_loc)/sig)) +\
                       (s>=top_loc)*1)
            # down-ramp
            else:
                npl = npl0*((s>top_loc)*(1/(1-2*(top_loc-s)/sig)) +\
                       (s<=top_loc)*1)
                
        # Xu PRL 2016 ramp shape 4
        elif shape.lower() == 'xu4':
            # define sigma parameter
            sig = hw/(np.sqrt(2)-1)
            # up-ramp
            if updn.lower() == 'up':
                npl = npl0*((s<top_loc)*(1/((1-(s-top_loc)/sig)**2)) +\
                       (s>=top_loc)*1)
            # down-ramp
            else:
                npl = npl0*((s>top_loc)*(1/((1-(top_loc-s)/sig)**2)) +\
                       (s<=top_loc)*1)
            
        # Xu PRL 2016 ramp shape 5
        elif shape.lower() == 'xu5':
            # define sigma parameter
            sig = hw/(2**(1/4)-1)
            # up-ramp
            if updn.lower() == 'up':
                npl = npl0*((s<top_loc)*(1/((1-(s-top_loc)/(2*sig))**4)) +\
                       (s>=top_loc)*1)
            # down-ramp
            else:
                npl = npl0*((s>top_loc)*(1/((1-(top_loc-s)/(2*sig))**4)) +\
                       (s<=top_loc)*1)


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

def insert_lens(plasma,lens_npl0,lens_L,lens_s0,add='yes'):
    """insert a thin plasma lens into plasma source"""
    s   = plasma["s"]
    npl = plasma["npl"]

    # find coarse start & stop indices
    istart  = np.argwhere(s>lens_s0)[0]-1
    istop   = np.argwhere(s>=lens_s0+lens_L)[0]
    
    if istart<0:
        istart=0
        istop=1
    
    icoarse = np.arange(istart,istop,dtype=np.int)
    
    # define fine indices, fine steps
    nfine   = int((s[istop]-s[istart])/(1e-6)) # m, new incremental size
    sfine   = np.linspace(s[istart],s[istop],nfine)
    ifstart = istart
    ifstop  = ifstart+nfine
    ifine   = np.arange(ifstart,ifstop,dtype=np.int)
    
    # remove old s values
    s = np.delete(s,icoarse)

    # insert new s values
    s = np.insert(s,istart,sfine)
    
    # create baseline npl values for fine steps
    slope   = (npl[istop]-npl[istart])/(s[istart+nfine]-s[istart])    
    nplfine = npl[istart] + slope*(s[ifine]-s[istart])
    
    # remove old npl values
    npl = np.delete(npl,icoarse)
    npl = np.insert(npl,istart,nplfine)

    # find lens start & stop indices
    ilstart = np.argwhere(s>=lens_s0)[0]
    ilstop  = np.argwhere(s>=lens_s0+lens_L)[0]
    ilens   = np.arange(ilstart,ilstop,dtype=np.int)
    
    # add lens density into npl
    if (add):
        npl[ilens] = npl[ilens]+lens_npl0
    else:
        npl[ilens] = lens_npl0
    
    # calculate new dgds array
    dgds0 = plasma["bulk"]["dgds0"]
    npl0  = plasma["bulk"]["npl0"]
    dgds  = calc_dgds(dgds0,npl0,npl)
    
    # modify plasma object
    plasma["s"]    = s
    plasma["npl"]  = npl
    plasma["dgds"] = dgds
    
    return plasma