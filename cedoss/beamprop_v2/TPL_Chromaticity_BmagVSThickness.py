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
len_arr = np.linspace(500, 7000, num) #original is 200 to 1200 for 1e18 cm-3
emit_arr = np.zeros(num)
betamin_arr = np.zeros(num)
tpl_f_arr = np.zeros(num)
centloc_arr = np.zeros(num)
betapromin_arr = np.zeros(num)
centbet_arr = np.zeros(num)

sigm_arr = np.zeros(num)
betam_arr = np.zeros(num)

position_error = 0

gammab = PProp.def_gamma
tpl_n = .5
betastar = .10 #0.00213065326633

delta = 0.01

for k in range(len(len_arr)):
    if k%2 == 0: print(k/len(len_arr)*100,"%");
    
    tpl_l = len_arr[k]
    tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
    tpl_f_arr[k] = tpl_f
    
    leftext = 1 #1
    rightext = 3 #3
    
    z_arr = np.linspace(-leftext*tpl_f, rightext*tpl_f, int((leftext+rightext)*tpl_f*1e6+1)*zmult) + (position_error / 1e6)
    n_arr = np.zeros(len(z_arr))+1e-9
    
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
    
    #emit_arr[k] = bmag_arr[-1]
    betapromin_arr[k] = min(betapro_arr)

    beam = PProp.GaussianBeam(beam_params, debug)
    dump = 10; cores = 4
    z_fine = np.copy(z_arr)*1e6
    PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)

    m = int(len(z_fine)/dump)
    emit_arr[k] = PProp.GetBmag(beam,m)
    sigm_arr[k] = PProp.GetSigmaMin(beam, m)
    betam_arr[k] = PProp.GetBetaMin(beam, m)


font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 12}

plt.rc('font', **font)
lwid = 2.0

kl = 1/tpl_f_arr
bw_arr = betastar * kl
dw = 0
sigmaE = 0.577 * delta
bmag_w2_arr = W2.CalcEmit(W2.ThinW2_Norm(bw_arr, dw),sigmaE)

#Thick
k = Foc.Calc_K(tpl_n*1e17, gammab)*100*100
bmag_w2_arr_thick = np.zeros(len(len_arr))
for x in range(len(len_arr)):
    w2 = W2.ThickW2_UnNormalized(k, len_arr[x]*1e-6, betastar, 0)
    bmag_w2_arr_thick[x] = W2.CalcEmit(w2, sigmaE)

projbeta_arr = np.zeros(len(len_arr))
projbeta_arr2 = np.zeros(len(len_arr))
for x in range(len(len_arr)):
    projbeta_arr[x] = W2.ProjBeta_UnNormalized(kl[x], betastar, 0, delta)
    b0 = betastar; a0 = 0; g0 = 1/betastar
    projbeta_arr2[x] = W2.ProjBetaCS_UnNormalized(k, len_arr[x]*1e-6, b0, a0, g0, delta)
    
projbetathick_arr = np.zeros(len(len_arr))
for x in range(len(len_arr)):
    projbetathick_arr[x] = W2.ProjBeta_Thick(k, len_arr[x]*1e-6, betastar, 0, delta)

plt.figure(figsize=(3.4*1.5, 2*1.5))
#plt.title("Emittance Growth vs Lens Thickness")
plt.plot(len_arr*1e-6*np.sqrt(k), emit_arr, 'r-', label = "Beam Propagation",linewidth = lwid)
plt.plot(len_arr*1e-6*np.sqrt(k), bmag_w2_arr, 'g-.', label = "Thin Calculated",linewidth = lwid)
plt.plot(len_arr*1e-6*np.sqrt(k), bmag_w2_arr_thick, 'b--', label = "Thick Calculated",linewidth = lwid)
plt.ylabel(r'$\epsilon_f/\epsilon_0$')
plt.xlabel(r'$\mathrm{Lens \ Thickness \ }\sqrt{K}l $')
plt.ylim(1.0,1.01)
plt.xlim(0.1, 1.5)
plt.grid(); plt.legend(); plt.show()

#plt.title("Minimum "+r'$\beta$'+" vs Lens Thickness")
#plt.plot(len_arr, betamin_arr*1e6, label="Measured")
#plt.plot(len_arr, Foc.Calc_BetaStar_DeltaOff(betastar, tpl_f_arr, 0)*1e6, label="Ideal Calculated")
#plt.plot(len_arr, projbeta_arr*1e6, label="Calculated "+r'$\beta_{pro}$')
#plt.plot(len_arr, betapromin_arr*1e6, label="Propagated "+r'$\beta_{pro}$')
#plt.ylabel(r'$\beta_{min}\mathrm{\ [\mu m]}$')
#plt.xlabel(r'$\mathrm{Lens \ Thickness \ [\mu m]}$')
#plt.grid(); plt.legend(); plt.show()

minloc = np.argmin(betamin_arr)
print("Minimum possible beta: ",betamin_arr[minloc]*1e6," um")
print("L = ",len_arr[minloc]," um")
print("B-mag = ",emit_arr[minloc])

#plt.plot(len_arr, np.sqrt(projbeta_arr*3e-6/gammab)*1e9, label="Thin Calculated "+r'$\sigma_r$')#, linewidth = lwid)
#plt.plot(len_arr, np.sqrt(projbetathick_arr*3e-6/gammab)*1e9, label="Thick Calculated "+r'$\sigma_r$')#, linewidth = lwid)
#plt.plot(len_arr, sigm_arr*1e9, label="Propagated "+r'$\sigma_r$')#, linewidth = lwid)
#plt.ylabel(r'$\sigma_r\mathrm{\ [nm]}$')
#plt.xlabel(r'$\mathrm{Lens \ Thickness \ [\mu m]}$')
#plt.grid(); plt.legend(); plt.show()

b0_arr = betastar
a0_arr = 0
p_arr = np.sqrt(1+np.square(b0_arr*sigmaE/tpl_f_arr))
sigp_arr = np.sqrt(b0_arr*3e-6/gammab)*p_arr/np.sqrt(np.square(p_arr)+np.square(a0_arr+b0_arr/tpl_f_arr))

##Nice plot of sigma r vs sqrt(k)l
plt.figure(figsize=(3.4*1.5, 2*1.5))
#plt.title("Beam spot size vs lens thickness")
plt.semilogy(len_arr*1e-6*np.sqrt(k), sigm_arr*1e9, 'r-', label="Beam Propagation", linewidth = lwid)
plt.semilogy(len_arr*1e-6*np.sqrt(k), np.sqrt(projbeta_arr*3e-6/gammab)*1e9, 'g-.', label="Thin Calculated", linewidth = lwid)
plt.semilogy(len_arr*1e-6*np.sqrt(k), np.sqrt(projbetathick_arr*3e-6/gammab)*1e9, 'b--', label="Thick Calculated", linewidth = lwid)
plt.semilogy(len_arr*1e-6*np.sqrt(k), np.sqrt(Foc.Calc_BetaStar_DeltaOff(betastar, tpl_f_arr, 0)*3e-6/gammab)*1e9, 'k:', label="Ideal", linewidth = lwid)
#plt.plot(len_arr*1e-6*np.sqrt(k), np.sqrt(Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(k,len_arr*1e-6,betastar, 0)*3e-6/gammab)*1e9, label="Thick Ideal "+r'$\sigma_r$')
#plt.semilogy(len_arr*1e-6*np.sqrt(k), sigp_arr*1e9, 'k--',label="Chen 1989 "+r'$\sigma_r$')
plt.ylabel(r'$\sigma_r\mathrm{\ [nm]}$')
plt.xlabel(r'$\mathrm{Lens \ Thickness \ }\sqrt{K}l $')
plt.grid(); plt.legend(); plt.show()

#plt.plot(len_arr, projbeta_arr, 'b-')
#plt.plot(len_arr, projbeta_arr2, 'r--')
#plt.show()
"""#Was testing new thick waist val and loc eqns
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
"""