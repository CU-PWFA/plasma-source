#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 14:01:20 2018

Functions to calculate analytic approximations for emittance growth

@author: chris
"""

import numpy as np

def ThinW2_Norm(bn, dn):
    return np.square(np.square(bn)+np.square(dn))/np.square(bn)

def ThinW2(kl, b, d):
    return np.square(kl)*(np.square(np.square(b)+np.square(d))/np.square(b))

def CalcEmit(w2, sigmaE):
    return 1 + 1/2*w2*np.square(sigmaE)

def ThickW2_UnNormalized(k, l, b, d):
    return ThickW2(l*np.sqrt(k), b*np.sqrt(k), d*np.sqrt(k))
    
#These are normalized to sqrt(k)
def ThickW2(l, b, d):
    first = (-d*l*np.cos(l)**2+(1+d*l)*np.sin(l)**2+.5*(d-l+b**2*l+d**2*l)*np.sin(2*l))**2 / ((b**2+d**2)*np.cos(l)**2+np.sin(l)**2+d*np.sin(2*l))**2
    second = .5*l*(1-b**2-d**2)*np.cos(l)**2-.5*l*(1-b**2-d**2)*np.sin(l)**2-.25*(1+b**2+d**2+4*d*l)*np.sin(2*l)
    thirdn = (-d*np.cos(l)**2+d*np.sin(l)**2-.5*(1-b**2-d**2)*np.sin(2*l))*(-d*l*np.cos(l)**2+(1+d*l)*np.sin(l)**2+.5*(d-l+b**2*l+d**2*l)*np.sin(2*l))
    thirdd = (b**2+d**2)*np.cos(l)**2+np.sin(l)**2+d*np.sin(2*l)
    return first + 1/b**2*(second - thirdn / thirdd)**2