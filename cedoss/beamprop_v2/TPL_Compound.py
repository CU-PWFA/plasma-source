#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 13:24:24 2018

Compound Lens

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

tpl_f1 = 0.05
tpl_f2 = 0.01
tpl_n = 10.
#tpl_l1 = 73.69 #um
#tpl_l2 = tpl_l1
#tpl_d = .0075 #m
tpl_l1 = Foc.Calc_Square_Lens(tpl_n*1e17, tpl_f1*100, gammab)
tpl_l2 = Foc.Calc_Square_Lens(tpl_n*1e17, tpl_f2*100, gammab)
tpl_d = 0.03 #m

tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l1, gammab)/100
tpl_f_plus = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l1, gammab*1.01)/100
tpl_f_mnus = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l1, gammab*0.99)/100
extend_domain = 1
z_arr = np.linspace(-1*tpl_f*extend_domain, 2*tpl_f*extend_domain, int(3*tpl_f*extend_domain*1e6+1)*zmult)
n_arr = np.zeros(len(z_arr))

dump = 10
cores = 4

betastar = .10 #0.00213065326633
waist_loc = 0.
tpl_offset = waist_loc
z_offset = -z_arr[0]
z_arr = z_arr + z_offset

tpl_1loc = tpl_offset-.5*tpl_d
tpl_2loc = tpl_offset+.5*tpl_d
tpl_1loc = tpl_offset
tpl_2loc = tpl_offset + tpl_d

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
argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_1loc*1e6, tpl_n, tpl_l1, debug)
argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_2loc*1e6, tpl_n, tpl_l2, debug)

fig, ax1 = plt.subplots()
plt.title("Beta evolution w/ "+r'$f_1=$'+'{:.3f}'.format(tpl_f1*100)+r'$\,\mathrm{cm}$'+", "+r'$f_2=$'+'{:.3f}'.format(tpl_f2*100)+r'$\,\mathrm{cm}$')
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
fig, ax3 = plt.subplots()
plt.title("Beta evolution w/ "+r'$f_1=$'+'{:.3f}'.format(tpl_f1*100)+r'$\,\mathrm{cm}$'+", "+r'$f_2=$'+'{:.3f}'.format(tpl_f2*100)+r'$\,\mathrm{cm}$')
ax3.set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'b')
ax3.tick_params('y', colors = 'b')
ax3.set_xlabel('z [cm]')

center0 = int((tpl_2loc+tpl_f)/(z_arr[-1]/(len(z_arr)+1)))
crange0 = 100*zmult
betacent0 = np.zeros(3)
centloc0 = np.zeros(3)
for i in range(len(e_spec)):
#Make beam and bulk plasma just as in single_pass
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                   beta_offset=waist_loc, plasma_start=z_offset,
                                                   gamma=gammab * e_spec[i])
    beam = PProp.GaussianBeam(beam_params, debug)
    
    beta = PProp.Calc_CSParams(beam_params, n_arr, z_arr)[0]
    betacent0[i] = beta[center0]
    ax3.plot((z_arr[center0-crange0:center0+crange0]-tpl_f)*1e2, np.array(beta[center0-crange0:center0+crange0])*1e2, colors[i], label=r'$\gamma/\gamma_{b}$' + " = "+str(e_spec[i]))
    
ax4 = ax3.twinx()
ax4.plot((z_arr[center0-crange0:center0+crange0]-tpl_f)*1e2, n_arr[center0-crange0:center0+crange0], 'k-')
ax4.set_ylabel(r'$n\,\mathrm{[10^{17}cm^{-3}]}$',color = 'k')
ax4.tick_params('y', colors = 'k')
ax3.grid(); ax3.legend(); plt.show()
print("beta at 2nd lens: ",betacent0)

###############################################################################

fig, ax5 = plt.subplots()
plt.title("Beta evolution w/ "+r'$f_1=$'+'{:.3f}'.format(tpl_f1*100)+r'$\,\mathrm{cm}$'+", "+r'$f_2=$'+'{:.3f}'.format(tpl_f2*100)+r'$\,\mathrm{cm}$')
ax5.set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'k')
ax5.tick_params('y', colors = 'k')
ax5.set_xlabel('z [cm]')

maxbetacomp = np.zeros(len(z_arr))
center = -1.
crange = 120*zmult
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
    for j in range(len(maxbetacomp)):
        if beta[j] > maxbetacomp[j]:
            maxbetacomp[j] = beta[j]
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

f_eff = 1/(1/tpl_f1+1/tpl_f2-tpl_d/tpl_f1/tpl_f2)
beta_pred = betastar/(1+(betastar/f_eff)**2)
z_v = tpl_d+(betastar/f_eff-tpl_d*betastar/f_eff/tpl_f1-tpl_d/betastar*(1-tpl_d/tpl_f2))/(
        betastar/f_eff**2+1/betastar*(1-tpl_d/tpl_f2)**2)
print("")
print("Effective Focal Length: ",f_eff)
print("Predicted Waist Value:  ",beta_pred)
print("     Error to sim [%]:  ",(betacent[0]-beta_pred)/beta_pred*100)
print("Predicted Waist Loc:    ",z_v)
print("   Error to sim [%]:    ",(centloc[0]-z_v)/z_v*100)
print("Minimum Max Beta: ",min(maxbetacomp))

###############################################################################

fig, ax6 = plt.subplots()
plt.title("Beta evolution w/ "+r'$f_1=$'+'{:.3f}'.format(tpl_f1*100)+r'$\,\mathrm{cm}$'+", "+r'$f_2=$'+'{:.3f}'.format(tpl_f2*100)+r'$\,\mathrm{cm}$')
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

beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar, beta_offset=waist_loc,
                                                plasma_start=z_offset, gamma=gammab)
gb_arr, beta_arr, alpha_arr, gamma_arr, bmag_arr = PProp.Calc_Proj_CSParams(beam_params, n_arr, z_arr, 0.01)

plt.title("B-mag vs z")
plt.ylabel("B-mag")
plt.xlabel("z [m]")
plt.plot(z_arr, bmag_arr)
plt.grid(); plt.show()
print("Final B-mag: ",bmag_arr[-1])