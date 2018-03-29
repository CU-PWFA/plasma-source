#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:39:45 2018

@author: chris
"""
import numpy as np

n_gam_fac = 56476779.72
cs = [0.329, 1.54, -0.0170]
cz = [-3.33, -4.47, 0.0897]
kb0 = 212.5

def CalcApprox(gamma_set, beta_set, kb_set):
    n_calc = n_gam_fac*kb_set**2*gamma_set
    sig_calc = cs[0]*beta_set**2*kb_set/kb0 + cs[1]*beta_set + cs[2]*kb0/kb_set
    waist_calc = cz[0]*beta_set**2*kb_set/kb0 + cz[1]*beta_set + cz[2]*kb0/kb_set
    return [gamma_set, beta_set, n_calc, sig_calc, waist_calc, 5*sig_calc]

def CalcApprox_Dens(gamma_set, beta_set, dens_set):
    kb = 1.33e-4*np.sqrt(dens_set / gamma_set)
    beta_n = beta_set * kb
    print("beta_n : ", beta_n)
    print("--Should be between 8.86 - 106--")
    
    sig_calc = (1.55e-3*beta_n**2 + 1.54*beta_n - 3.61)/kb
    waist_calc = (-1.57e-2*beta_n**2 - 4.47*beta_n + 19.1)/kb
    return [gamma_set, beta_set, dens_set, sig_calc, waist_calc, 5*sig_calc]

if __name__ == '__main__':
    gamma_set = 19569.5
    beta_set = 0.10
    dens_set = 5e16
#    kb_set = 212.5
    
    #vals = CalcApprox(gamma_set,beta_set,kb_set)
    vals = CalcApprox_Dens(gamma_set, beta_set, dens_set)
    
    print("gamma_set = " +str(vals[0]))
    print("beta_set = " + str(vals[1]))
    print("dens = " + str(vals[2]))
    print("hwup_set = " + str(vals[3]))
    print("waist_set = " + str(vals[4]))
    print("L_up_set = " +str(vals[5]))