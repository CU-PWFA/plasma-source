#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 12:33:24 2018

Uses macro-particles to numerically caluclate emittance growth and beam rms size
slower but more authentic.

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
zmult=10

num = 101
d_arr = np.linspace(-0.30, -0.00, num)
emit_arr = np.zeros(num)
betamin_arr = np.zeros(num)
betapromin_arr= np.zeros(num)
sigm_arr = np.zeros(num)
betam_arr = np.zeros(num)
gammaf_arr = np.zeros(num)

gammab = PProp.def_gamma
tpl_n = 10.
tpl_l = 800
betastar = .10 #0.00213065326633

delta = 0.01

for k in range(len(d_arr)):
    if k%10 == 0: print(k/len(d_arr)*100,"%");
    position_error = d_arr[k] * 1e6 #um
    
    tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
    
    leftext = 1 #1
    rightext = 3 #3
    
    z_arr = np.linspace(-leftext*tpl_f, rightext*tpl_f, int((leftext+rightext)*tpl_f*1e6+1)*zmult) + (position_error / 1e6)
    n_arr = np.zeros(len(z_arr))#+1e-9
    
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
    
    gb_arr, beta_arr, alpha_arr, gamma_arr, bmag_arr, betapro_arr = PProp.Calc_Proj_CSParams(beam_params, n_arr, z_arr, delta)
    
    betapromin_arr[k] = min(betapro_arr)
    emit_arr[k] = bmag_arr[-1]
    
    beam = PProp.GaussianBeam(beam_params, debug)
    dump = 10; cores = 4
    z_fine = np.copy(z_arr)*1e6
    PProp.PropBeamPlasma(beam, argon, z_fine, dump, cores, debug)

    m = int(len(z_fine)/dump)
    emit_arr[k] = PProp.GetBmag(beam,m)
    sigm_arr[k] = PProp.GetSigmaMin(beam, m)
    betam_arr[k] = PProp.GetBetaMin(beam, m)
    #print(sigm_arr[k])
    #print(emit_arr[k])
    gammaf_arr[k] = PProp.GetFinalGamma(beam,m)

#tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
#gammaf=np.average(gammaf_arr)
#ave_factor = 1+(0.5*(gammaf-gammab)/gammab)
#tpl_f=tpl_f*ave_factor

kl = 1/tpl_f
bw = betastar * kl
dw_arr = d_arr * kl
sigmaE = 0.577 * delta # sqrt(1/3) = 0.577
w2_arr = W2.ThinW2_Norm(bw, dw_arr)
bmag_w2_arr = W2.CalcEmit(w2_arr,sigmaE)

#Thick
l = tpl_l*1e-6
k = 1/(l*tpl_f)
bmag_w2_arr_thick = np.zeros(len(d_arr))
for x in range(len(d_arr)):
    w2 = W2.ThickW2_UnNormalized(k, l, betastar, d_arr[x])
    bmag_w2_arr_thick[x] = W2.CalcEmit(w2, sigmaE)

projbeta_arr = np.zeros(len(d_arr))
projbeta_arr2 = np.zeros(len(d_arr))
for x in range(len(d_arr)):
    projbeta_arr[x] = W2.ProjBeta_UnNormalized(kl, betastar, d_arr[x], delta)
    b0 = betastar+d_arr[x]**2/betastar;  a0 = -d_arr[x]/betastar;  g0 = 1/betastar
    projbeta_arr2[x] = W2.ProjBetaCS_UnNormalized(k, l, b0, a0, g0, delta)
    
projbetathick_arr = np.zeros(len(d_arr))
for x in range(len(d_arr)):
    #projbetathick_arr[x] = W2.ProjBeta_Thick(k, l, betastar, d_arr[x], delta)
    projbetathick_arr[x] = W2.ProjBeta_Thick_Gauss(k, l, betastar, d_arr[x], sigmaE)

plt.title("B-mag vs Lens-Waist Separation for L = "+str(tpl_l)+r'$\ \mu m$')
plt.plot(d_arr*1e2, emit_arr, 'r-', label = "Beam Propagation")
plt.plot(d_arr*1e2, bmag_w2_arr, 'g--', label = "Analytic Thin")
plt.plot(d_arr*1e2, bmag_w2_arr_thick, 'b--', label = "Analytic Thick")
plt.ylabel("B-mag")
plt.xlabel("Lens-Waist Separation d [cm]")
plt.grid(); plt.legend(); plt.show()

plt.title("Minimum "+r'$\beta$'+" vs Lens-Waist Separation for L = "+str(tpl_l)+r'$\ \mu m$')
#plt.plot(d_arr*1e2, betamin_arr*1e6, label="Measured")
plt.plot(d_arr*1e2, betam_arr*1e6, label="Statistical")
plt.plot(d_arr*1e2, Foc.Calc_BetaStar_DeltaOff(betastar, tpl_f, d_arr)*1e6, label="Ideal Calculated")
plt.plot(d_arr*1e2, projbeta_arr*1e6, label="Calculated "+r'$\beta_{pro}$')
plt.plot(d_arr*1e2, projbetathick_arr*1e6, label="Thick Calculated "+r'$\beta_{pro}$')
plt.plot(d_arr*1e2, betapromin_arr*1e6, label="Propagated "+r'$\beta_{pro}$')
plt.ylabel(r'$\beta_{min}\mathrm{\ [\mu m]}$')
plt.xlabel("Lens-Waist Separation d [cm]")
#plt.ylim([.1*min(betamin_arr*1e6), 1*max(betamin_arr*1e6)])
plt.grid(); plt.legend(); plt.show()

minloc = np.argmin(betamin_arr)
print("Minimum possible beta: ",betamin_arr[minloc]*1e6," um")
print("d = ",d_arr[minloc]*100," cm")
print("B-mag = ",emit_arr[minloc])

b0_arr = betastar + np.square(d_arr)/betastar
a0_arr = -d_arr/betastar
p_arr = np.sqrt(1+np.square(b0_arr*sigmaE/tpl_f))
sigp_arr = np.sqrt(b0_arr*3e-6/gammab)*p_arr/np.sqrt(np.square(p_arr)+np.square(a0_arr+b0_arr/tpl_f))

##Nice plot of sigma r vs sqrt(k)d
plt.title("Beam spot size vs lens-waist separation")
#plt.plot(d_arr*np.sqrt(k), np.sqrt(projbeta_arr*3e-6/gammab)*1e9, label="Thin Calculated "+r'$\sigma_r$')
plt.plot(d_arr*np.sqrt(k), np.sqrt(projbetathick_arr*3e-6/gammab)*1e9, 'b--', label="Thick Calculated "+r'$\sigma_r$')
plt.plot(d_arr*np.sqrt(k), sigm_arr*1e9, 'r-', label="Propagated "+r'$\sigma_r$')
plt.plot(d_arr*np.sqrt(k), np.sqrt(Foc.Calc_BetaStar_DeltaOff(betastar, tpl_f, d_arr)*3e-6/gammab)*1e9, 'g-', label="Ideal "+r'$\sigma_r$')
#plt.plot(d_arr*np.sqrt(k), np.sqrt(Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(k,tpl_l*1e-6,betastar, d_arr)*3e-6/gammab)*1e9, label="Thick Ideal "+r'$\sigma_r$')
#plt.plot(d_arr*np.sqrt(k), sigp_arr*1e9, 'k--',label="Chen 1989 "+r'$\sigma_r$')
plt.ylabel(r'$\sigma_r\mathrm{\ [nm]}$')
plt.xlabel(r'$\mathrm{Lens-Waist \ Separation \ }\sqrt{K}d $')
plt.grid(); plt.legend(); plt.show()

#This was to confirm new betapro formula worked
#plt.plot(d_arr, projbeta_arr, 'b-')
#plt.plot(d_arr, projbeta_arr2, 'r--')
#plt.show()

"""
#Below is production figure

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}

plt.rc('font', **font)
lwid = 3.0
#fig, ax1 = plt.subplots()
#plt.plot(d_arr*1e2, np.sqrt(projbeta_arr*3e-6*bmag_w2_arr/gammab)*1e9, 'b-', label="Thin Approx.", linewidth = lwid)
plt.plot(d_arr*1e2, np.sqrt(projbeta_arr*3e-6/gammab)*1e9, 'b-', label="Thin Approx.", linewidth = lwid)
#plt.plot(d_arr*1e2, np.sqrt(projbetathick_arr*3e-6*bmag_w2_arr_thick/gammab)*1e9, 'c-', label="Thick Regime", linewidth = lwid)
plt.plot(d_arr*1e2, np.sqrt(projbetathick_arr*3e-6/gammab)*1e9, 'c-', label="Thick Regime", linewidth = lwid)
plt.plot(d_arr*1e2, sigm_arr*1e9, 'r--', label="Numerical", linewidth = lwid)
plt.plot(d_arr*1e2, np.sqrt(Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(k,tpl_l*1e-6,betastar, d_arr)*3e-6/gammab)*1e9, 'g-', label="Ideal "+r'$\sigma_E=0$', linewidth = lwid)
plt.ylabel(r'$\sigma_r\mathrm{\ [nm]}$')
plt.xlabel(r'$\mathrm{Lens-Waist \ Separation} \ d \ \mathrm{[cm]}$')
plt.ylim([60, 125])
#fig.set_size_inches(6,5)
plt.grid(); plt.legend(loc=(.22,.63)); plt.show()
"""