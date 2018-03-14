#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:39:45 2018

@author: chris
"""
n_gam_fac = 56476779.72

def Getkb(kb_option):
    if kb_option == 1:
        return 94.0
    if kb_option == 2:
        return 212.5
    if kb_option == 3:
        return 297.3
    if kb_option == 4:
        return 210.2

def CalcApprox(gamma_set, beta_set, kb_option=2):
    if kb_option==1: 
        kb_set = 94.0
        n_calc = n_gam_fac*kb_set**2*gamma_set
        sig_calc = 1.48*beta_set - 0.0209
        waist_calc = -3.42*beta_set**2 - 3.30*beta_set + 0.0678
        
    elif kb_option==2: 
        kb_set = 212.48975158758344
        n_calc = n_gam_fac*kb_set**2*gamma_set
        sig_calc = 1.65*beta_set - 0.0246
        waist_calc = -3.33*beta_set**2 - 4.47*beta_set + 0.0897
        
    elif kb_option==3: 
        kb_set = 297.254100056
        n_calc = n_gam_fac*kb_set**2*gamma_set
        sig_calc = 1.71*beta_set - 0.0192
        waist_calc = -3.92*beta_set**2 - 4.67*beta_set + 0.0670
    elif kb_option==4: 
        kb_set = 210.190389885
        n_calc = n_gam_fac*kb_set**2*gamma_set
        sig_calc = 1.63*beta_set - 0.0172
        waist_calc = -2.68*beta_set**2 - 4.57*beta_set + 0.0761
    return [gamma_set, beta_set, n_calc, sig_calc, waist_calc, 5*sig_calc]

if __name__ == '__main__':
    gamma_set = 19569.5
    beta_set = 0.10
    kb_opt = 2
    
    vals = CalcApprox(gamma_set,beta_set,kb_opt)
    
    print("gamma_set = " +str(vals[0]))
    print("beta_set = " + str(vals[1]))
    print("dens = " + str(vals[2]))
    print("hwup_set = " + str(vals[3]))
    print("waist_set = " + str(vals[4]))
    print("L_up_set = " +str(vals[5]))