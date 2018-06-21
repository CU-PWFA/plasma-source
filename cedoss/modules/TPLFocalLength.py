#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 12:39:09 2017

Calculates the focal length as an integral of density

Here focal refers to the 'focal length' = 1/KL
The actual waist position is at z=focal(1-beta_f/beta_i) for simple case,
and is given by Calc_Waist_Delta for offset cases

@author: chris
"""
import numpy as np
import scipy.integrate as Int

gam_def = 19569.5
const_def = 5.648e15

#These are CGS and used in recent TPL analysis:
def Calc_K(den, gam):
    e = 4.8032e-10
    m = 9.1094e-28
    c = 2.9979e10
    K = (2*np.pi*den*e**2)/(gam*m*c**2)
    return K

def Calc_Focus_KLength(K, plen):
    return 1/(K*plen)

#densities in cm^-3
#focal legnths in cm
#lens thicknesses in um
#Integrates over a density profile to calculate the focal length (n cm)
def Calc_Focus(den,axis,gam=gam_def,const=const_def,scale=1e17):
    focal = gam*const/(Int.simps(den*scale,axis))
    print(str(focal) + " cm focus")
    return focal

#Given density and length, what is the focal length (in cm)
    #Had a bunch of issues, so either SI or cm & um
def Calc_Focus_Square_CM_UM(den, plen, gam):
    focal = 5.648e15*gam/den/plen
    return focal

def Calc_Focus_Square_SI(den, plen, gam):
    focal = 5.648e7*gam/den/plen
    return focal

#Given focal length and density, how long is TPL (in um)
def Calc_Square_Lens(den, focal, gam=gam_def, const=const_def):
    plen = gam*const/den/focal
    return plen
    
#Given initial and final beta, what is required f (in cm)
def Calc_Target_Focal(beta_i, beta_f):
    focal = np.sqrt(beta_i**2 * beta_f/(beta_i - beta_f))
    return focal

#Given focal length and initial beta, what is final beta (in cm)
def Calc_BetaStar(beta_i, focal):
    beta_f = beta_i/(1+beta_i**2/focal**2)
    return beta_f

#Given focal length and initial beta and offset delta, what is final beta
def Calc_BetaStar_DeltaOff(beta_i, focal, delta=0):
    beta_f = beta_i/((beta_i/focal)**2+(1-delta/focal)**2)
    return beta_f

#Given focla length, initial beta, and offset delta; what is distance from lens to waist
def Calc_WaistPos_DeltaOff(beta_i, focal, delta=0):
    z = ((1/focal)*(beta_i**2 + delta**2) - delta)/((beta_i/focal)**2+(1-delta/focal)**2)
    return z

if __name__ == '__main__':
    beta_i = .10
    beta_f = 0.011693
    n0=5e16*5
    gam_set = gam_def
    
    focal = Calc_Target_Focal(beta_i,beta_f)
    #focal = 0.015
    f_len = Calc_Square_Lens(n0, focal*100, gam=gam_set)
    print(f_len," : lens length [um]")
    print(focal*100," : focal length[cm]")
    print(750*np.sqrt(gam_set/n0)*1e6," : max length [um]")
    
    #print(Calc_Focus_Square_CM_UM(1e18, 74.91, 19565.5)/100)
    #print(Calc_Focus_Square_SI(5e16, 100e-6, 2e4))