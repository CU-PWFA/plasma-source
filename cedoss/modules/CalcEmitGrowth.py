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

def ProjBeta_UnNormalized(kl, b, d, delta):
    return ProjBeta(b*kl, d*kl, delta)/kl
    
#Normalized to k*l.  Valid for thin lenses but it works well for thick too.
def ProjBeta(b, d, delta):
    first = 2*b*delta*(b**2+(d-1)**2)**2
    second = 2*delta/(delta**2-1)*(b**6*(delta**2-2)+(delta**2-2)*(d-1)**2*d**4+b**2*(d-1)*(3*d-1)*(delta**2-1+(delta**2-2)*d**2)+b**4*(-3+4*d-6*d**2+delta**2*(3+d*(3*d-2))))
    third = -4*(b**2+(d-1)*d)*(b**2-d*(b**2+d**2)+(b**2+d**2)**2)*np.arctanh(delta)
    return (second + third) / first

#This approach used an unnecesary assumption
"""
def ProjBeta_UnNormalized_OLD(kl, b, d, delta):
    return ProjBeta_OLD(b*kl, d*kl, delta)/kl
    
#Normalized to k*l.  Valid for thin lenses.
def ProjBeta_OLD(b, d, delta):
    first = 2*delta*b*(b**2 + (b**2 + d**2)**2)
    second = (b**4 + 3*b**2*d +d**3 - d**4)*(np.arctan(b/(d-1+delta))-np.arctan(b/(d-1-delta)))
    third = (b**3*(d-1) + b*d**3)*np.log((b**2+(1-d-delta)**2)/(b**2+(1-d+delta)**2))
    return 1/(2*delta*b**2)*(first + (b**2+d**2-d)*(second - third))
"""