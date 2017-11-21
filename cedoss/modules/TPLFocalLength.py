#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 12:39:09 2017

Calculates the focal length as an integral of density

@author: chris
"""
import numpy as np
import scipy.integrate as Int

gam_def = 2e4
const_def = 5.648e15

def Calc_Focus(den,axis,gam=gam_def,const=const_def,scale=1e17):
    focal = gam*const/(Int.simps(den*scale,axis))
    print(str(focal) + " cm focus")
    return focal