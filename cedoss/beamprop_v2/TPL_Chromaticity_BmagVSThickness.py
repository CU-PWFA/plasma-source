#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 14:08:03 2018

Copy of TPL_Chromaticity_LargeOffset built specifically to measure emittance
growth vs lens thickness

@author: chris
"""

import sys
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import CalcEmitGrowth as W2

path = '/home/chris/Desktop/BeamProp/Placeholder'
debug = 0
zmult=1

num = 101
len_arr = np.linspace(800, 1200, num)
emit_arr = np.zeros(num)
betamin_arr = np.zeros(num)
tpl_f_arr = np.zeros(num)
centloc_arr = np.zeros(num)
betapromin_arr = np.zeros(num)
centbet_arr = np.zeros(num)
position_error = 0

gammab = PProp.def_gamma
tpl_n = 10.
betastar = .10 #0.00213065326633

delta = 0.01

for k in range(len(len_arr)):
    if k%10 == 0: print(k/len(len_arr)*100,"%");
    
    tpl_l = len_arr[k]
    tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
    tpl_f_arr[k] = tpl_f
    
    leftext = 1 #1
    rightext = 3 #3
    
    z_arr = np.linspace(-leftext*tpl_f, rightext*tpl_f, int((leftext+rightext)*tpl_f*1e7+1)*zmult) + (position_error / 1e6)
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
    argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6 + position_error + tpl_l/2, tpl_n, tpl_l, debug)
    
    ###############################################################################
    
    maxbetacomp = np.zeros(len(z_arr))
    center = -1.
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
    
    betamin_arr[k] = min(maxbetacomp)
    centloc_arr[k] = centloc[0]-tpl_l*1e-6
    centbet_arr[k] = betacent[0]
    ###############################################################################
    
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar, beta_offset=waist_loc,
                                                    plasma_start=z_offset, gamma=gammab)
    
    gb_arr, beta_arr, alpha_arr, gamma_arr, bmag_arr, betapro_arr  = PProp.Calc_Proj_CSParams(beam_params, n_arr, z_arr, delta)
    
    emit_arr[k] = bmag_arr[-1]
    betapromin_arr[k] = min(betapro_arr)

kl = 1/tpl_f_arr
bw_arr = betastar * kl
dw = 0
sigmaE = 0.57 * delta
bmag_w2_arr = W2.CalcEmit(W2.ThinW2_Norm(bw_arr, dw),sigmaE)

#Thick
k = Foc.Calc_K(tpl_n*1e17, gammab)*100*100
bmag_w2_arr_thick = np.zeros(len(len_arr))
for x in range(len(len_arr)):
    w2 = W2.ThickW2_UnNormalized(k, len_arr[x]*1e-6, betastar, 0)
    bmag_w2_arr_thick[x] = W2.CalcEmit(w2, sigmaE)

projbeta_arr = np.zeros(len(len_arr))
for x in range(len(len_arr)):
    projbeta_arr[x] = W2.ProjBetaAlt_UnNormalized(kl[x], betastar, 0, delta)

plt.title("B-mag vs Lens Thickness")
plt.plot(len_arr, emit_arr, label = "Beam Propagation")
plt.plot(len_arr, bmag_w2_arr, label = "Analytic Thin")
plt.plot(len_arr, bmag_w2_arr_thick, label = "Analytic Thick")
plt.ylabel("B-mag")
plt.xlabel(r'$\mathrm{Lens \ Thickness \ [\mu m]}$')
plt.grid(); plt.legend(); plt.show()

plt.title("Minimum "+r'$\beta$'+" vs Lens Thickness")
plt.plot(len_arr, betamin_arr*1e6, label="Measured")
plt.plot(len_arr, Foc.Calc_BetaStar_DeltaOff(betastar, tpl_f_arr, 0)*1e6, label="Ideal Calculated")
plt.plot(len_arr, projbeta_arr*1e6, label="Calculated "+r'$\beta_{pro}$')
plt.plot(len_arr, betapromin_arr*1e6, label="Propagated "+r'$\beta_{pro}$')
plt.ylabel(r'$\beta_{min}\mathrm{\ [\mu m]}$')
plt.xlabel(r'$\mathrm{Lens \ Thickness \ [\mu m]}$')
plt.grid(); plt.legend(); plt.show()

minloc = np.argmin(betamin_arr)
print("Minimum possible beta: ",betamin_arr[minloc]*1e6," um")
print("L = ",len_arr[minloc]," um")
print("B-mag = ",emit_arr[minloc])

plt.plot(len_arr, np.sqrt(projbeta_arr*3e-6*bmag_w2_arr/gammab)*1e9, label="Calculated "+r'$\sigma_r$')
plt.plot(len_arr, np.sqrt(projbeta_arr*3e-6*bmag_w2_arr_thick/gammab)*1e9, label="Thick Calculated "+r'$\sigma_r$')
plt.plot(len_arr, np.sqrt(betapromin_arr*3e-6*emit_arr/gammab)*1e9, label="Propagated "+r'$\sigma_r$')
plt.ylabel(r'$\sigma_r\mathrm{\ [nm]}$')
plt.xlabel("Lens-Waist Separation d [cm]")
plt.grid(); plt.legend(); plt.show()

k = Foc.Calc_K(tpl_n*1e17, gammab)*100*100

plt.title("Thick Waist Location vs Simulations")
plt.plot(len_arr, centloc_arr, label = "Center Loc.")
plt.plot(len_arr, tpl_f_arr, label = "Focal Length")
plt.plot(len_arr, Foc.Calc_ThickWaistPos_DeltaOff_UnNormalized(k, len_arr*1e-6, betastar, 0), label = "Thick Calculated")
plt.grid(); plt.legend(); plt.show()

plt.title("Thick Waist Value vs Simulations")
plt.plot(len_arr, centbet_arr, label = "Waist Value")
plt.plot(len_arr, Foc.Calc_BetaStar_DeltaOff(betastar, tpl_f_arr, 0), label = "Thin Eqn")
plt.plot(len_arr, Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(k, len_arr*1e-6, betastar, 0), label = "Thick Eqn")
plt.grid(); plt.legend(); plt.show()