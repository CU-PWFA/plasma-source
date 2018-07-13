#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 13:04:09 2018

Calculate's mathematica's solution to lens-waist separation matching

@author: chris
"""

import numpy as np
import sys
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

bi = .061608 #initial beta
bf = .061608 #final beta / matching condition
zm = 0.10 #distance from inital waist to final waist - set by facility and mathcing condition

if bi == bf:
    print("Using bi = bf case:")
    kl = 4*zm/(4*bi**2 + zm**2)
    d = .5*zm
    
    print("KL = ",kl,"m^-1")
    print("d  = ",d,"m")
    
    gamma = Foc.gam_def
    tpl_n = 0.5
    
    tpl_l = Foc.Calc_Square_Lens(tpl_n*1e17, 1/kl*100, gamma)
    print("L  = ",tpl_l,"um")
    
else:
    """#old formulas
    kl = (2*bf*bi*zm + np.sqrt(bf*bi*(bf+bi)**2*(bi**2-2*bf*bi+bf**2+zm**2))) / (bf*bi*(bi**2+2*bf*bi+bf**2+zm**2))
    d = (bf*bi*zm + bi**2*zm - np.sqrt(bf*bi*(bf+bi)**2*(bi**2-2*bf*bi+bf**2+zm**2))) / (bi**2-bf**2)
    """
    kl = (2*bf*bi*zm + (bi+bf)*np.sqrt(bf*bi*((bi-bf)**2+zm**2))) / (bf*bi*((bi+bf)**2+zm**2))
    
    d =  (bi*zm-np.sqrt(bi*bf*((bi-bf)**2+zm**2)))/(bi-bf)
    
    print("zm = ",zm, "m")
    print("KL = ",kl,"m^-1")
    print("d  = ",d,"m")
    
    gamma = Foc.gam_def
    tpl_n = 0.5
    
    tpl_l = Foc.Calc_Square_Lens(tpl_n*1e17, 1/kl*100, gamma)
    print("L  = ",tpl_l,"um")
    print()
    
    if bf > bi:
        x = 1.30
        print("Lower limit for zm when x = ",x)
        zmin = np.sqrt(bi*(x*bf-bi))+bf*np.sqrt(x-1)
        print("zmin = ",zmin)
    
    #It appears taking the other sign allos for negative KL - unphysical
    """
    kl = (2*bf*bi*zm - np.sqrt(bf*bi*(bf+bi)**2*(bi**2-2*bf*bi+bf**2+zm**2))) / (bf*bi*(bi**2+2*bf*bi+bf**2+zm**2))
    
    d = (-1*bf*bi*zm -bi**2*zm - np.sqrt(bf*bi*(bf+bi)**2*(bi**2-2*bf*bi+bf**2+zm**2))) / (bf**2-bi**2)
    
    print()
    print("***Taking other sign:")
    print()
    
    print("KL",kl)
    print("d ",d)
    """