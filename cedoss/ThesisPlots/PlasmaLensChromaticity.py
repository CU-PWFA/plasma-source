#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 16:58:36 2023

Chromatic analysis of TPL

This version is a copy of TPL_Chromaticity for Thesis fig 2.2

@author: chris
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from beamprop_v2 import BeamPropFuncs as PProp

path = '/home/chris/Desktop/BeamProp/GasCellTest'
debug = 0

zmult=1

gammab = PProp.def_gamma
tpl_n = 10.

tpl_f = 0.01#0.10
tpl_l = Foc.Calc_Square_Lens(tpl_n*1e17, tpl_f*100, gammab)

#tpl_l = 1400

tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
tpl_f_plus = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab*1.01)/100
tpl_f_mnus = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab*0.99)/100

leftext = 1 #1
rightext = 3 #3

z_arr = np.linspace(-leftext*tpl_f, rightext*tpl_f, int((leftext+rightext)*tpl_f*1e6+1)*zmult)
n_arr = np.zeros(len(z_arr))

dump = 10
cores = 4

betastar = .10 #0.00213065326633
waist_loc = 0.
tpl_offset = waist_loc
z_offset = -z_arr[0]
z_arr = z_arr + z_offset

position_error = 0#-7*tpl_f*1e6 #um

e_spec = np.array([0, -0.01, 0.01]) + 1.0
colors = np.array(['g-','r-','b-'])
arrlist = np.array([])

beam_params0 = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                   beta_offset=waist_loc, plasma_start=z_offset,
                                                   gamma=gammab)
argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                                   nset = tpl_n)
argon_params['Z'] = z_arr[-1]*1e6
argon_params['Nz']= len(z_arr)
argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6 + position_error, tpl_n, tpl_l, debug)

fig, ax1 = plt.subplots(figsize=(11.5,3.5),nrows=1, ncols=2)

#plt.title("Beta function evolution at "+r'$f=$'+'{:.3f}'.format(tpl_f*100)+r'$\,\mathrm{cm}$')
ax1[0].set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'b')
ax1[0].tick_params('y', colors = 'b')
ax1[0].set_xlabel('z [cm]')
ax1[0].set_ylim([-0.05,20.05])
beta0 = PProp.Calc_CSParams(beam_params0, np.zeros(len(z_arr)), z_arr)[0]
ax1[0].plot((z_arr-tpl_f)*1e2, np.array(beta0)*1e2, 'b--',label=r'$\beta_{vac}$')
ax1[0].text(-.93,18.5,"(a)")

for i in range(len(e_spec)):
#Make beam and bulk plasma just as in single_pass
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                   beta_offset=waist_loc, plasma_start=z_offset,
                                                   gamma=gammab * e_spec[i])
    beta = PProp.Calc_CSParams(beam_params, n_arr, z_arr)[0]
    ax1[0].plot((z_arr-tpl_f)*1e2, np.array(beta)*1e2, colors[i], label=r'$\gamma/\gamma_{b}$' + " = "+str(e_spec[i]))
    
ax2 = ax1[0].twinx()
ax2.plot((z_arr-tpl_f)*1e2, n_arr, 'k-')
ax2.set_ylabel(r'$n\,\mathrm{[10^{17}cm^{-3}]}$',color = 'k')
ax2.tick_params('y', colors = 'k')
ax2.set_ylim([0,12.1])
ax1[0].grid(); ax1[0].legend(loc=4); #plt.show()

###############################################################################

#fig, ax5 = plt.subplots()
#plt.subplot(1,2,2)
#plt.title("Beta function evolution at "+r'$f=$'+'{:.3f}'.format(tpl_f*100)+r'$\,\mathrm{cm}$')
ax1[1].set_ylabel(r'$\beta\,\mathrm{[mm]}$', color = 'k')
ax1[1].tick_params('y', colors = 'k')
ax1[1].set_xlabel('z [cm]')

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
    ax1[1].plot((z_arr[center-crange:center+crange]-tpl_f)*1e2, np.array(beta)*1e2*10, colors[i], label=r'$\gamma/\gamma_{b}$' + " = "+str(e_spec[i]))

ax1[1].axvline(x=(centloc[0])*100, c='g', ls='--')
ax1[1].axvline(x=(centloc[1])*100, c='r', ls='--')
ax1[1].axvline(x=(centloc[2])*100, c='b', ls='--')
ax1[1].text(.971,1.087,"(b)")

ax1[1].grid(); ax1[1].legend(loc=1); 
plt.tight_layout()
plt.show()