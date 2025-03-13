#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 12:24:56 2019

For the lithium oven plasma source, single beam propagations to demonstrate
the ideal matching case and some TPL solutions.

@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import sys
sys.path.insert(0, "../../python")
from lens import profile

debug = 1
path = '/home/chris/Desktop/BeamProp/testGaussian'

tanhlens_option = 0

case = 2
if case == 1:#Ideal Matching
    tpl_l = 0
    tpl_offset = 0
    tpl_n = 0 
    betastar = 0.04214 #m
    waistloc = 0.0357 #m
    
tpl_foc_dist = 0.6684
if case == 2:#10cm with lens far upstream
    tpl_l = 55.1273022719
    tpl_n = 0.1
    d = -1.08197312858
    tpl_offset = -tpl_foc_dist
    
    betastar = 0.1
    waistloc = tpl_offset - d
if case == 3:#10cm with lens far downstream
    tpl_l = 257.700528841
    tpl_n = 0.1
    d = 1.08197312858
    tpl_offset = -tpl_foc_dist
    
    betastar = 0.1
    waistloc = tpl_offset - d

zw = 0.8168
tpl_foc_dist = zw - 0.0357
if case == 11:#10cm with lens far upstream - post E300 Conference
    tpl_l = 52.787916222855195
    tpl_n = 0.09
    d = -1.2559530948862472
    tpl_offset = -tpl_foc_dist
    
    betastar = 0.1
    waistloc = tpl_offset - d
if case == 12:#10cm with lens far upstream - post E300 Conference
    tpl_l = 107.05292551328377
    tpl_n = 0.09
    d = -2.77256034790439
    tpl_offset = -tpl_foc_dist
    
    betastar = 0.5
    waistloc = tpl_offset - d

zw = 0.6009
tpl_foc_dist = zw - 0.0357
if case == 21:#10cm with lens far upstream - post E300 Conference
    tpl_l = 19.005416805373816
    tpl_n = 0.33999
    d = -0.9225367278532891
    tpl_offset = -tpl_foc_dist
    
    betastar = 0.1
    waistloc = tpl_offset - d
if case == 22:#10cm with lens far upstream - post E300 Conference
    tpl_l = 38.62821117036434
    tpl_n = 0.33999
    d = -2.0137949426868835
    tpl_offset = -tpl_foc_dist
    
    betastar = 0.5
    waistloc = tpl_offset - d

zw = 0.6009
tpl_foc_dist = zw - 0.0357 
if case == 31: #With 5Torr Helium buffer gas, only need a 3.87 um lens given by these tanh params
    tpl_a = 14.522400931076138
    tpl_b = 11.900068868957758
    tpl_n0 = 0.1334178042556154 * 1.67
    d = -0.922536727853289
    tpl_offset = -tpl_foc_dist
    tanhlens_option = 1
    betastar = 0.1
    waistloc = tpl_offset - d

argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = 0., scaledown = 1)

leftext = 0.7811#1 #m
ovenext = 0.5
dz = .1 #um
zarr = np.arange((-ovenext-leftext)*1e6, ovenext*1e6+dz, dz)
argon_params['Nz'] = len(zarr)
argon_params['Z'] = zarr[-1]-zarr[0]

lioven = profile.lithium_oven_profile(z=zarr, center=0, ne0=0.34)
start = lioven[0]/1e6
nfunc = lioven[2]
narr=nfunc(zarr)

if debug == 1:
    plt.plot(zarr,narr)
    plt.grid(); plt.show()
    
relset = ovenext + leftext + start
waistloc = waistloc + relset
tpl_offset = tpl_offset + relset

#print(tpl_n0)

if tanhlens_option == 0:
    argon = PProp.CustomPlasma_ThinPlasmaLens(argon_params, narr, tpl_offset*1e6, tpl_n, tpl_l, debug)
else:
    argon = PProp.CustomPlasma_ThinPlasmaLens_Tanh(argon_params, narr, tpl_offset*1e6, [tpl_a,tpl_b,tpl_n0], debug)
n = argon.nez
zarr = zarr - zarr[0]
z = zarr/1e6

#Make beam and bulk plasma just as in single_pass
beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waistloc, plasma_start=0)
beam_params_matched = PProp.ReturnDefaultElectronParams(path, beta_star=0.04214,
                                                       beta_offset=0.0357+relset, plasma_start=0)
endindex = np.argmax(n)
bmag = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
        #bmag = PProp.Calc_Bmag(beam_params,n, z)
print("bmag",bmag)
if case == 1:
    betaf,alphaf,gammaf,gbf = PProp.Plot_CSEvo(beam_params, n, z, 0, legend_loc = 0)
else:
    betaf,alphaf,gammaf,gbf = PProp.Plot_CSEvo_MatchedCompare(beam_params, beam_params_matched,n, z, 0, legend_loc = 0)