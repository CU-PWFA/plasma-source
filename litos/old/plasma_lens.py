#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  4 14:50:40 2017

@author: litos
"""

import numpy as np

def plasma_lens(npl,dgds,s,args,dgds0=[0]):
    
    [np0,s0,Lp] = args
    
    n_lens = 0
    i_lo, = np.where(s<=s0)
    if (i_lo.any()):
        i_lo  = i_lo[-1]
        i_hi, = np.where(s>=s[i_lo]+Lp)
        if (i_hi.any()):
            i_hi  = i_hi[0]    
            n_lens = i_hi-i_lo
    
#    print(args)
#    print(i_lo)
#    print(i_hi)
#    print(n_lens)
    if n_lens>0:
        npl[i_lo:i_hi]  = np0*np.ones(n_lens)
        dgds[i_lo:i_hi] = dgds0*np.ones(n_lens)
    
    return [npl, dgds]