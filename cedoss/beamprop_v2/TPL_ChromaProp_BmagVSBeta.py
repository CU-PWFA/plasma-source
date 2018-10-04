#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 10:25:21 2018

Copy of TPL_ChromaProp_LargeOffset built specifically to measure emittance
growth vs initial beta (no lens-waist separation)

@author: chris
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 14:08:03 2018



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
"""
num = 101
beta_arr = np.linspace(0.05, 0.50, num)
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
tpl_n = 10.
tpl_l = 800

delta = 0.01

for k in range(len(beta_arr)):
    if k%10 == 0: print(k/len(beta_arr)*100,"%");
    
    betastar = beta_arr[k]
    tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
    tpl_f_arr[k] = tpl_f
    
    leftext = 1 #1
    rightext = 2 #3
    
    z_arr = np.linspace(-leftext*tpl_f, rightext*tpl_f, int((leftext+rightext)*tpl_f*1e6+1)*zmult) + (position_error / 1e6)
    n_arr = np.zeros(len(z_arr))+1e-9
    
    waist_loc = 0.
    tpl_offset = waist_loc
    z_offset = -z_arr[0]
    z_arr = z_arr + z_offset
    
    e_spec = np.array([0, -delta, delta]) + 1.0
    
    argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                                       nset = tpl_n)
    argon_params['Z'] = z_arr[-1]*1e6
    argon_params['Nz']= len(z_arr)
    argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
    argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6 + position_error + tpl_l/2, tpl_n, tpl_l, debug)
    
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar, beta_offset=waist_loc,
                                                    plasma_start=z_offset, gamma=gammab)
    
    gb_arr, beta_arr_differentone, alpha_arr, gamma_arr, bmag_arr, betapro_arr  = PProp.Calc_Proj_CSParams(beam_params, n_arr, z_arr, delta)
    
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


#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 15}
#
#plt.rc('font', **font)
lwid = 3.0
"""
kl = 1/tpl_f
bw_arr = beta_arr * kl
dw = 0
sigmaE = 0.577 * delta
bmag_w2_arr = W2.CalcEmit(W2.ThinW2_Norm(bw_arr, dw),sigmaE)

#Thick
k = Foc.Calc_K(tpl_n*1e17, gammab)*100*100
bmag_w2_arr_thick = np.zeros(len(beta_arr))
for x in range(len(beta_arr)):
    w2 = W2.ThickW2_UnNormalized(k, tpl_l*1e-6, beta_arr[x], 0)
    bmag_w2_arr_thick[x] = W2.CalcEmit(w2, sigmaE)

projbeta_arr = np.zeros(len(beta_arr))
projbeta_arr2 = np.zeros(len(beta_arr))
for x in range(len(beta_arr)):
    projbeta_arr[x] = W2.ProjBeta_UnNormalized(kl, beta_arr[x], 0, delta)
    b0 = beta_arr[x]; a0 = 0; g0 = 1/beta_arr[x]
    projbeta_arr2[x] = W2.ProjBetaCS_UnNormalized(k, tpl_l*1e-6, b0, a0, g0, delta)
    
projbetathick_arr = np.zeros(len(beta_arr))
for x in range(len(beta_arr)):
    projbetathick_arr[x] = W2.ProjBeta_Thick(k, tpl_l*1e-6, beta_arr[x], 0, delta)

plt.title("Emittance Growth vs Lens Thickness")
plt.plot(beta_arr*np.sqrt(k), bmag_w2_arr, label = "Thin Calculated")#,linewidth = lwid)
plt.plot(beta_arr*np.sqrt(k), bmag_w2_arr_thick, label = "Thick Calculated")#,linewidth = lwid)
plt.plot(beta_arr*np.sqrt(k), emit_arr, label = "Beam Propagation")#,linewidth = lwid)
plt.ylabel(r'$\epsilon_f/\epsilon_0$')
plt.xlabel(r'$\mathrm{Initial \ Beta \ }\sqrt{K}\beta_i^* $')
plt.grid(); plt.legend(); plt.show()

plt.plot(beta_arr, np.sqrt(projbeta_arr*3e-6/gammab)*1e9, label="Thin Calculated "+r'$\sigma_r$')#, linewidth = lwid)
plt.plot(beta_arr, np.sqrt(projbetathick_arr*3e-6/gammab)*1e9, label="Thick Calculated "+r'$\sigma_r$')#, linewidth = lwid)
plt.plot(beta_arr, sigm_arr*1e9, label="Propagated "+r'$\sigma_r$')#, linewidth = lwid)
plt.ylabel(r'$\sigma_r\mathrm{\ [nm]}$')
plt.xlabel(r'$\mathrm{Initial \ Beta \ }\beta_i^* \ \mathrm{[m]}$')
plt.grid(); plt.legend(); plt.show()

b0_arr = beta_arr
a0_arr = 0
p_arr = np.sqrt(1+np.square(b0_arr*sigmaE/tpl_f_arr))
sigp_arr = np.sqrt(b0_arr*3e-6/gammab)*p_arr/np.sqrt(np.square(p_arr)+np.square(a0_arr+b0_arr/tpl_f_arr))

##Nice plot of sigma r vs sqrt(k)l
plt.title("Beam spot size vs initial beta")
#plt.semilogy(beta_arr*np.sqrt(k), np.sqrt(projbeta_arr*3e-6/gammab)*1e9, label="Thin Calculated "+r'$\sigma_r$')
plt.semilogy(beta_arr, np.sqrt(projbetathick_arr*3e-6/gammab)*1e9, label="Thick Calculated "+r'$\sigma_r$')#, linewidth = lwid)
plt.semilogy(beta_arr, sigm_arr*1e9, label="Propagated "+r'$\sigma_r$')
plt.semilogy(beta_arr, np.sqrt(Foc.Calc_BetaStar_DeltaOff(beta_arr, tpl_f, 0)*3e-6/gammab)*1e9, label="Thin Ideal "+r'$\sigma_r$')
#plt.plot(len_arr*1e-6*np.sqrt(k), np.sqrt(Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(k,len_arr*1e-6,betastar, 0)*3e-6/gammab)*1e9, label="Thick Ideal "+r'$\sigma_r$')
#plt.semilogy(beta_arr*np.sqrt(k), sigp_arr*1e9, 'k--',label="Chen 1989 "+r'$\sigma_r$')
plt.ylabel(r'$\sigma_r\mathrm{\ [nm]}$')
plt.xlabel(r'$\mathrm{Initial \ Beta \ }\beta_i^* \ \mathrm{[m]}$')
plt.grid(); plt.legend(); plt.show()

plt.plot(beta_arr, projbeta_arr, 'b-')
plt.plot(beta_arr, projbeta_arr2, 'r--')
plt.show()