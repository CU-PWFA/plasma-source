#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 12:35:49 2018

Calculates values from Oide paper

@author: chris
"""

import numpy as np
import scipy.integrate as Int

pi = np.pi
"""#SI values
h = 6.6261e-34
m = 9.1094e-31
c = 2.9979e8
comp = h/(m*c)
re = 2.8179e-15
"""
h = 6.6261e-27
m = 9.1094e-28
c = 2.9979e10
comp = 3.8616e-11
re = 2.8179e-13

def Get_KLls_WRONNNGNG(lens_width, focal_length):
    K = 1/focal_length/lens_width
    L = lens_width
    ls = focal_length - 1/2*lens_width
    return [K,L,ls]

def Get_ls(lens_width, focal_length):
    return focal_length - 1/2*lens_width

def Inner_F(lens, focs, x):
    return Int.quad(lambda y: np.square(np.sin(y)+focs*np.cos(y)), 0, x)[0]

#Where KLls = [K,L,ls] and is returned by Get_KLls
def F_Oide(KLls):
    lens = np.sqrt(KLls[0])*KLls[1]
    focs = np.sqrt(KLls[0])*KLls[2]
    return Int.quad(lambda x: np.power(np.abs(np.sin(x)+focs*np.cos(x)),3)*
                      np.square(Inner_F(lens,focs, x)), 0, lens)[0]

def Calc_SigMin(F_val, emit):
    return (7/5)**(1/2) * (275/(3*np.sqrt(6*pi))*re*comp*F_val)**(1/7) * emit**(5/7)

def Calc_SigOide(F_val, emit, gam, beta_f):
    return np.sqrt(beta_f*emit/gam + 110/(3*np.sqrt(6*pi))*re*comp*gam**5*F_val*(emit/gam/beta_f)**(5/2))

def Calc_SigIdeal(F_val, emit, gam, beta_f):
    return np.sqrt(beta_f*emit/gam)

def Calc_BetaMin(F_val, emit, gam):
    return (275/(3*np.sqrt(6*pi))*re*comp*F_val)**(2/7) * gam * emit**(3/7)