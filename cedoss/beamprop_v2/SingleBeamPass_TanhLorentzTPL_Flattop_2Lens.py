#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 11:38:24 2019

@author: chris
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 15:39:05 2019

TPL matching into hard cuttof flattop plasma source

This version uses a PTPL with DoubleTanh-Lorentzian density
distribution, and will cater towards the PWFA example in my
PTPL Paper 1.

Specifically, we want both an entrance and an exit ramp to focus to 10cm beta

@author: chris
"""

import sys
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

debug = 1
path = '/home/chris/Desktop/BeamProp/testGaussian'
gamma = PProp.def_gamma

sighw = 0.0
tpl_n = 0.5
tpl_l = 994.99
zvac = -0.02084 - .5*tpl_l*1e-6
betastar = .10
    
a = 495.627005532
b = 60.5967508142
n_0 = 0.5

A = 3.780e3*1000000
gamma = 4812.9*1000000
x_0 = 0

fit_tanh = [a,b,n_0]
fit_lorentz = [A, gamma, x_0]

case = 1

if case == 1:
    tpl_l2 = 1115.8
    zvac2 = 0.10 + 0.02159 + .5*tpl_l2*1e-6
    #zvac2 = 0.10 + 0.02091 + .5*tpl_l2*1e-6
    #zvac2 = 0.10 + 0.02070 + .5*tpl_l2*1e-6
    
    a2 = 555.991426567
    b2 = 53.9334556975
    fit_tanh2 = [a2,b2,n_0]

if case == 2:
    tpl_l2 = 1142.85
    #zvac2 = 0.10 + 0.02111 + .5*tpl_l2*1e-6
    zvac2 = 0.10 + 0.02171 + .5*tpl_l2*1e-6
    
    a2 = 569.514448709
    b2 = 52.6388627684
    fit_tanh2 = [a2,b2,n_0]
    
if case == 3:
    tpl_l2 = tpl_l
    zvac2 = 0.10 + 0.02246 + .5*tpl_l2*1e-6
    
    a2 = a
    b2 = b
    fit_tanh2 = [a2,b2,n_0]

z0 = 0.1

argon_params = PProp.ReturnDefaultPlasmaParams(path, sigma_hw = sighw, plasma_start = z0, scaledown = 10)
argon = PProp.GaussianRampPlasma(argon_params, debug)

n = argon.nez
z = np.linspace(0,argon_params['Z'],len(n))/1e6

dump = 10
cores = 4

ramp_end = argon_params['z0']/1e6
endindex = np.nonzero(z>ramp_end)[0][0]

#Make beam and bulk plasma just as in single_pass
beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                               beta_offset=zvac, plasma_start=z0)
beam = PProp.GaussianBeam(beam_params, debug)
argon = PProp.CustomPlasma_ThinPlasmaLens_TanhLorentz(argon_params, n, zvac*1e6, fit_tanh, fit_lorentz, debug)
argon = PProp.CustomPlasma_ThinPlasmaLens_TanhLorentz(argon_params, n, zvac2*1e6, fit_tanh2, fit_lorentz, debug)

z_fine = np.copy(z)*1e6
PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)

m = int(len(z_fine)/dump)

PProp.PlotEmittance(beam,z_fine,m)
PProp.PlotGamma(beam,z_fine,m)

print("Bmag BP: ",PProp.GetBmag(beam,m))
bmagc = PProp.Calc_Bmag(beam_params,n[:endindex], z[:endindex])
print("Bmag CS: ",bmagc)

betaf, alphaf, gammaf, gbf = PProp.Plot_CSEvo(beam_params, n, z, z0, legend_loc = 10)

print(1/gammaf[-1])
#set_sta = 8000; set_fin = 12000
#PProp.Plot_CSEvo(beam_params, n, z, z0, legend_loc = 0, subset = [set_sta, set_fin])