#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 15:12:05 2018

Copy of TPL_Chromaticity_LargeOffset built specifically to measure emittance
growth vs lens-waist separation

@author: chris
"""

import sys
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

path = '/home/chris/Desktop/BeamProp/Placeholder'
debug = 0
zmult=1

num = 101
d_arr = np.linspace(-0.20, -0.00, num)
emit_arr = np.zeros(num)
betamin_arr = np.zeros(num)

gammab = PProp.def_gamma
tpl_n = 10.
tpl_l = 1105
betastar = .10 #0.00213065326633

delta = 0.01

for k in range(len(d_arr)):
    if k%10 == 0: print(k/len(d_arr)*100,"%");
    position_error = d_arr[k] * 1e6 #um
    
    tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
    
    leftext = 1 #1
    rightext = 3 #3
    
    z_arr = np.linspace(-leftext*tpl_f, rightext*tpl_f, int((leftext+rightext)*tpl_f*1e6+1)*zmult) + (position_error / 1e6)
    n_arr = np.zeros(len(z_arr))
    
    waist_loc = 0.
    tpl_offset = waist_loc
    z_offset = -z_arr[0]
    z_arr = z_arr + z_offset
    
    e_spec = np.array([0, -delta, delta]) + 1.0
    
    beam_params0 = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=z_offset,
                                                       gamma=gammab)
    argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                                       nset = tpl_n)
    argon_params['Z'] = z_arr[-1]*1e6
    argon_params['Nz']= len(z_arr)
    argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
    argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6 + position_error, tpl_n, tpl_l, debug)
    
    ###############################################################################
    
    maxbetacomp = np.zeros(len(z_arr))
    center = -1.
    crange = 200*zmult
    betacent = np.zeros(3)
    centloc = np.zeros(3)
    for i in range(len(e_spec)):
    #Make beam and bulk plasma just as in single_pass
        beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=z_offset,
                                                       gamma=gammab * e_spec[i])
        beta = PProp.Calc_CSParams(beam_params, n_arr, z_arr)[0]
        if center < 0:
            center = np.argmin(beta)
        for j in range(len(maxbetacomp)):
            if beta[j] > maxbetacomp[j]:
                maxbetacomp[j] = beta[j]
        centloc[i] = z_arr[np.argmin(beta)]-tpl_f
        betacent[i] = beta[center]
        beta = beta[center-crange:center+crange]
    
    betamin_arr[k] = min(maxbetacomp)
    
    ###############################################################################
    
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar, beta_offset=waist_loc,
                                                    plasma_start=z_offset, gamma=gammab)
    
    gb_arr, beta_arr, alpha_arr, gamma_arr, bmag_arr = PProp.Calc_Proj_CSParams(beam_params, n_arr, z_arr, delta)
    
    emit_arr[k] = bmag_arr[-1]

kl = 1/tpl_f
bw = betastar * kl
dw_arr = d_arr * kl
sigmaE = 0.6 * delta
bmag_w2_arr = 1 + 1/2*(np.square(bw**2+np.square(dw_arr))/bw**2)*sigmaE**2

plt.title("B-mag vs Lens-Waist Separation for L = "+str(tpl_l)+r'$\ \mu m$')
plt.plot(d_arr*1e2, emit_arr, label = "Beam Propagation")
plt.plot(d_arr*1e2, bmag_w2_arr, label = "Analytical Approximation")
plt.ylabel("B-mag")
plt.xlabel("Lens-Waist Separation d [cm]")
plt.grid(); plt.show()

plt.title("Minimum "+r'$\beta$'+" vs Lens-Waist Separation for L = "+str(tpl_l)+r'$\ \mu m$')
plt.plot(d_arr*1e2, betamin_arr*1e6)
plt.ylabel(r'$\beta_{min}\mathrm{\ [\mu m]}$')
plt.xlabel("Lens-Waist Separation d [cm]")
plt.grid(); plt.show()

minloc = np.argmin(betamin_arr)
print("Minimum possible beta: ",betamin_arr[minloc]*1e6," um")
print("d = ",d_arr[minloc]*100," cm")
print("B-mag = ",emit_arr[minloc])