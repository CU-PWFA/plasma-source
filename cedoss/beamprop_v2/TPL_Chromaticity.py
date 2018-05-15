#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 04 10:52:25 2018

Chromatic analysis of TPL

@author: chris
"""

import sys
import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

path = '/home/chris/Desktop/BeamProp/GasCellTest'
debug = 0

zmult=1

gammab = PProp.def_gamma
tpl_n = 10.
tpl_l = 184.21#74.91
tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
tpl_f_plus = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab*1.01)/100
tpl_f_mnus = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab*0.99)/100

z_arr = np.linspace(-1*tpl_f, 2*tpl_f, int(3*tpl_f*1e6+1)*zmult)
n_arr = np.zeros(len(z_arr))

dump = 10
cores = 4

betastar = .10 #0.00213065326633
waist_loc = 0.
tpl_offset = waist_loc
z_offset = -z_arr[0]
z_arr = z_arr + z_offset

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
argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6, tpl_n, tpl_l, debug)

fig, ax1 = plt.subplots()
plt.title("Beta function evolution at "+r'$f=$'+'{:.3f}'.format(tpl_f*100)+r'$\,\mathrm{cm}$')
ax1.set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'b')
ax1.tick_params('y', colors = 'b')
ax1.set_xlabel('z [cm]')
ax1.set_ylim([-0.05,20.05])
beta0 = PProp.Calc_CSParams(beam_params0, np.zeros(len(z_arr)), z_arr)[0]
ax1.plot((z_arr-tpl_f)*1e2, np.array(beta0)*1e2, 'b--',label=r'$\beta_{vac}$')

for i in range(len(e_spec)):
#Make beam and bulk plasma just as in single_pass
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                   beta_offset=waist_loc, plasma_start=z_offset,
                                                   gamma=gammab * e_spec[i])
    beam = PProp.GaussianBeam(beam_params, debug)
    

    
    beta = PProp.Calc_CSParams(beam_params, n_arr, z_arr)[0]
    ax1.plot((z_arr-tpl_f)*1e2, np.array(beta)*1e2, colors[i], label=r'$\gamma/\gamma_{b}$' + " = "+str(e_spec[i]))
    
ax2 = ax1.twinx()
ax2.plot((z_arr-tpl_f)*1e2, n_arr, 'k-')
ax2.set_ylabel(r'$n\,\mathrm{[10^{17}cm^{-3}]}$',color = 'k')
ax2.tick_params('y', colors = 'k')
ax1.grid(); ax1.legend(); plt.show()

###############################################################################
"""
fig, ax3 = plt.subplots()
plt.title("Beta function evolution at "+r'$f=$'+'{:.3f}'.format(tpl_f*100)+r'$\,\mathrm{cm}$')
ax3.set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'b')
ax3.tick_params('y', colors = 'b')
ax3.set_xlabel('z ['+r'$\mu$'+'m]')
#ax3.set_ylim([-0.05,20.05])
beta0 = PProp.Calc_CSParams(beam_params0, np.zeros(len(z_arr)), z_arr)[0]
center0 = np.argmin(beta0)
crange0 = 10*zmult
ax3.plot((z_arr[center0-crange0:center0+crange0]-tpl_f)*1e6, np.array(beta0[center0-crange0:center0+crange0])*1e2, 'b--',label=r'$\beta_{vac}$')

for i in range(len(e_spec)):
#Make beam and bulk plasma just as in single_pass
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                   beta_offset=waist_loc, plasma_start=z_offset,
                                                   gamma=gammab * e_spec[i])
    beam = PProp.GaussianBeam(beam_params, debug)
    
    beta = PProp.Calc_CSParams(beam_params, n_arr, z_arr)[0]
    ax3.plot((z_arr[center0-crange0:center0+crange0]-tpl_f)*1e6, np.array(beta[center0-crange0:center0+crange0])*1e2, colors[i], label=r'$\gamma/\gamma_{b}$' + " = "+str(e_spec[i]))
    
ax4 = ax3.twinx()
ax4.plot((z_arr[center0-crange0:center0+crange0]-tpl_f)*1e6, n_arr[center0-crange0:center0+crange0], 'k-')
ax4.set_ylabel(r'$n\,\mathrm{[10^{17}cm^{-3}]}$',color = 'k')
ax4.tick_params('y', colors = 'k')
ax3.grid(); ax3.legend(); plt.show()
"""
###############################################################################

fig, ax5 = plt.subplots()
plt.title("Beta function evolution at "+r'$f=$'+'{:.3f}'.format(tpl_f*100)+r'$\,\mathrm{cm}$')
ax5.set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'k')
ax5.tick_params('y', colors = 'k')
ax5.set_xlabel('z [cm]')

center = -1.
crange = 50*zmult
betacent = np.zeros(3)
centloc = np.zeros(3)
for i in range(len(e_spec)):
#Make beam and bulk plasma just as in single_pass
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                   beta_offset=waist_loc, plasma_start=z_offset,
                                                   gamma=gammab * e_spec[i])
    beam = PProp.GaussianBeam(beam_params, debug)
    
    beta = PProp.Calc_CSParams(beam_params, n_arr, z_arr)[0]
    if center < 0:
        center = np.argmin(beta)
    
    centloc[i] = z_arr[np.argmin(beta)]-tpl_f
    betacent[i] = beta[center]
    beta = beta[center-crange:center+crange]
    ax5.plot((z_arr[center-crange:center+crange]-tpl_f)*1e2, np.array(beta)*1e2, colors[i], label=r'$\gamma/\gamma_{b}$' + " = "+str(e_spec[i]))
"""
ax5.axvline(x=(centloc[0])*100, c='g')
ax5.axvline(x=(centloc[1])*100, c='r')
ax5.axvline(x=(centloc[2])*100, c='b')
"""
ax5.grid(); ax5.legend(); plt.show()
print("beta at center's waist: ",betacent)
print("waist location for beta: ",centloc)
print((z_arr[center]-tpl_f)*100)