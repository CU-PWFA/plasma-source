#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 14:01:20 2018

Functions to calculate analytic approximations for emittance growth

@author: chris
"""

import numpy as np
import scipy.integrate as Int


def ThinW2_Norm(bn, dn):
    return np.square(np.square(bn)+np.square(dn))/np.square(bn)

def ThinW2(kl, b, d):
    return np.square(kl)*(np.square(np.square(b)+np.square(d))/np.square(b))

def CalcEmit(w2, sigmaE):
    return np.sqrt(1+w2*np.square(sigmaE))

def CalcEmit_OLD(w2, sigmaE):
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

def ProjBetaCS_UnNormalized(k, l, b0, a0, g0, delta):
    return ProjBetaCS(l*np.sqrt(k), b0*np.sqrt(k), a0, g0/np.sqrt(k), delta)/np.sqrt(k)

def ProjBetaCS(l, b, a, g, delta):
    first = 2*delta*(g+l*(2*a+b*l))**2
    second = 2*delta/(delta**2-1)*((g+4*a*l)*(delta**2-1)+b*l**2*(a**2*(1-2*delta**2)+3*b*g*(-1+delta**2))+b**2*l**3*(2*a+b*l)*(-2+delta**2))
    third = -4*l*(a+b*l)*(1+a*b*l+b**2*l**2)*np.arctanh(delta)
    #print(first,second,third)
    #print((second + third) / first)
    return (second + third) / first

def ProjBeta_Thick(k, l, b, d, delta):
    rk = np.sqrt(k)
    rkl = rk*l
    srkl = np.sin(rkl)
    crkl = np.cos(rkl)
    
    numterm = -d*rk*crkl**2 - crkl*srkl + b**2*k*crkl*srkl + d**2*k*crkl*srkl + d*rk*srkl**2
    denterm = rk*(crkl**2-2*d*rk*crkl*srkl+b**2*k*srkl**2+d**2*k*srkl**2)
    
    t11 = numterm / denterm
    t21 = -t11*np.sqrt(k)
    t12 = 1/rk
    
    return Int.quad(lambda x: 1/b*(t11*np.cos(rkl/np.sqrt(x))+t12*np.sin(rkl/np.sqrt(x))/(1/np.sqrt(x)))**2+
                    (b+d**2/b)*(np.cos(rkl/np.sqrt(x))+t21*(1/np.sqrt(x))*np.sin(rkl/np.sqrt(x)))**2 + 
                    2*d/b * (t11*np.cos(rkl/np.sqrt(x))+t12*np.sin(rkl/np.sqrt(x))/(1/np.sqrt(x))) * 
                    (np.cos(rkl/np.sqrt(x))+t21*(1/np.sqrt(x))*np.sin(rkl/np.sqrt(x))), 1-delta, 1+delta)[0]/(2*delta)

def ProjBeta_Thick_Gauss(k, l, b, d, sigE):
    rk = np.sqrt(k)
    rkl = rk*l
    srkl = np.sin(rkl)
    crkl = np.cos(rkl)
    
    numterm = -d*rk*crkl**2 - crkl*srkl + b**2*k*crkl*srkl + d**2*k*crkl*srkl + d*rk*srkl**2
    denterm = rk*(crkl**2-2*d*rk*crkl*srkl+b**2*k*srkl**2+d**2*k*srkl**2)
    
    t11 = numterm / denterm
    t21 = -t11*np.sqrt(k)
    t12 = 1/rk
    
    return Int.quad(lambda x: (1/b*(t11*np.cos(rkl/np.sqrt(1+x))+t12*np.sin(rkl/np.sqrt(1+x))/(1/np.sqrt(1+x)))**2+
                    (b+d**2/b)*(np.cos(rkl/np.sqrt(1+x))+t21*(1/np.sqrt(1+x))*np.sin(rkl/np.sqrt(1+x)))**2 + 
                    2*d/b * (t11*np.cos(rkl/np.sqrt(1+x))+t12*np.sin(rkl/np.sqrt(1+x))/(1/np.sqrt(1+x))) * 
                    (np.cos(rkl/np.sqrt(1+x))+t21*(1/np.sqrt(1+x))*np.sin(rkl/np.sqrt(1+x))))*np.exp(-np.square(x)/(2*np.square(sigE)))
                    /(np.sqrt(2*np.pi)*sigE), -10*sigE, 10*sigE)[0]

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