#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 13:21:40 2018

Version of TPL_Chromaticity meant to handle large offsets

As it turned out, all I really needed to do was shift z_arr also by the amount 
'position_error', this way z_arr is already shifted but the rest of the code
understands that it is really offset so the beam is initialized differently

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

#tpl_f = 0.001
#tpl_l = Foc.Calc_Square_Lens(tpl_n*1e17, tpl_f*100, gammab)
tpl_l = 1500

position_error = -0.0 * 1e6#-7*tpl_f*1e6 #um

delta = 0.01

tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
tpl_f_plus = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab*1.01)/100
tpl_f_mnus = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab*0.99)/100

leftext = 1 #1
rightext = 3 #3

z_arr = np.linspace(-leftext*tpl_f, rightext*tpl_f, int((leftext+rightext)*tpl_f*1e6+1)*zmult) + (position_error / 1e6)
n_arr = np.zeros(len(z_arr))

dump = 10
cores = 4

betastar = .10 #0.00213065326633
waist_loc = 0.
tpl_offset = waist_loc
z_offset = -z_arr[0]
z_arr = z_arr + z_offset

e_spec = np.array([0, -delta, delta]) + 1.0
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
argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6 + position_error + tpl_l/2, tpl_n, tpl_l, debug)

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
    beta = PProp.Calc_CSParams(beam_params, n_arr, z_arr)[0]
    ax1.plot((z_arr-tpl_f)*1e2, np.array(beta)*1e2, colors[i], label=r'$\gamma/\gamma_{b}$' + " = "+str(e_spec[i]))
    
ax2 = ax1.twinx()
ax2.plot((z_arr-tpl_f)*1e2, n_arr, 'k-')
ax2.set_ylabel(r'$n\,\mathrm{[10^{17}cm^{-3}]}$',color = 'k')
ax2.tick_params('y', colors = 'k')
ax1.grid(); ax1.legend(); plt.show()

###############################################################################

fig, ax5 = plt.subplots()
plt.title("Beta function evolution at "+r'$f=$'+'{:.3f}'.format(tpl_f*100)+r'$\,\mathrm{cm}$')
ax5.set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'k')
ax5.tick_params('y', colors = 'k')
ax5.set_xlabel('z [cm]')

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
    ax5.plot((z_arr[center-crange:center+crange]-tpl_f)*1e2, np.array(beta)*1e2, colors[i], label=r'$\gamma/\gamma_{b}$' + " = "+str(e_spec[i]))

ax5.axvline(x=(centloc[0])*100, c='g')
ax5.axvline(x=(centloc[1])*100, c='r')
ax5.axvline(x=(centloc[2])*100, c='b')

ax5.grid(); ax5.legend(); plt.show()

print("beta at center's waist: ",betacent)
print("waist location for beta: ",centloc)
print((z_arr[center]-tpl_f)*100)

print("Minimum Max Beta: ",min(maxbetacomp))

###############################################################################

fig, ax6 = plt.subplots()
plt.title("Beta function evolution at "+r'$f=$'+'{:.3f}'.format(tpl_f*100)+r'$\,\mathrm{cm}$')
ax6.set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'k')
ax6.tick_params('y', colors = 'k')
ax6.set_xlabel('z [cm]')

crangef = 120*zmult
#centerf = len(z_arr)-2*crangef
centerf = -1
betaf = np.zeros(3)
dbetaf = np.zeros(3)
for i in range(len(e_spec)):
#Make beam and bulk plasma just as in single_pass
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                   beta_offset=waist_loc, plasma_start=z_offset,
                                                   gamma=gammab * e_spec[i])
    beam = PProp.GaussianBeam(beam_params, debug)
    beta = PProp.Calc_CSParams(beam_params, n_arr, z_arr)[0]
    if centerf < 0:
        centerf = np.where(np.array(beta)[int(.5*len(beta)):]>0.11)[0][0]+int(.5*len(beta))
    betaf[i]=beta[centerf]
    dbetaf[i]=(beta[centerf+1]-beta[centerf-1])/(z_arr[centerf+1]-z_arr[centerf-1])

    beta = beta[centerf-crangef:centerf+crangef]
    ax6.plot((z_arr[centerf-crangef:centerf+crangef]-tpl_f)*1e2, np.array(beta)*1e2, colors[i], label=r'$\gamma/\gamma_{b}$' + " = "+str(e_spec[i]))

ax6.grid(); ax6.legend(); plt.show()
print("Beta final vals: ",betaf)
print("Beta spread [%]: ",100*(betaf[1]-betaf[2])/betaf[0])
print("Beta derivative: ",dbetaf)
print("Der. spread [%]: ",100*(dbetaf[1]-dbetaf[2])/dbetaf[0])

###############################################################################

print("Attempting with projected")
beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar, beta_offset=waist_loc,
                                                plasma_start=z_offset, gamma=gammab)

gb_arr, beta_arr, alpha_arr, gamma_arr, bmag_arr = PProp.Calc_Proj_CSParams(beam_params, n_arr, z_arr, delta)

plt.title("B-mag vs z")
plt.ylabel("B-mag")
plt.xlabel("z [m]")
plt.plot(z_arr, bmag_arr)
plt.grid(); plt.show()
print("Final B-mag: ",bmag_arr[-1])
"""
beam = PProp.GaussianBeam(beam_params, debug)
argon = PProp.CustomPlasma(argon_params, n_arr+.00000000001, 1)
PProp.PropBeamPlasma(beam, argon, z_arr*1e6, dump, cores, debug)
m = int(len(z_arr)/dump)-1
PProp.PlotEmittance(beam,z_arr*1e6,m)
print("Bmag: ",PProp.GetBmag(beam,m))
"""