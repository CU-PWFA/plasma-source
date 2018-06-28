#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 14:34:15 2018

Loop over thickness - find transition from thin

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

num = 201
thick_arr = np.linspace(50,4000,num)
betas_arr = np.zeros(num)
waist_arr = np.zeros(num)
focus_arr = np.zeros(num)
focus2_arr = np.zeros(num)
betap_arr = np.zeros(num)

for i in range(len(thick_arr)):
    zmult=1
    gammab = PProp.def_gamma
    tpl_n = 1.
    tpl_l = thick_arr[i]
    tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
    
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
    
    argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                                       nset = tpl_n)
    argon_params['Z'] = z_arr[-1]*1e6
    argon_params['Nz']= len(z_arr)
    argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
    argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6 + position_error, tpl_n, tpl_l, debug)
    
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                       beta_offset=waist_loc, plasma_start=z_offset,
                                                       gamma=gammab)
    beta = PProp.Calc_CSParams(beam_params, n_arr, z_arr)[0]
    center = np.argmin(beta)
    centloc = z_arr[center]-tpl_f
    betacent = beta[center]
    
    waist_arr[i] = centloc
    betas_arr[i] = betacent
    focus_arr[i] = tpl_f
    focus2_arr[i] = tpl_f*(1-betacent/betastar)
    betap_arr[i] = Foc.Calc_BetaStar(betastar, tpl_f)

plt.title("Waist Locations:  Propagation vs Theory")
plt.semilogy(thick_arr, waist_arr*100, label = "Waist locations")
#plt.semilogy(thick_arr, focus_arr*100)
plt.semilogy(thick_arr, focus2_arr*100, label = "f*(1-b_f/b_i)")
plt.xlabel(r'$\mathrm{Lens \ Thickness \ [\mu m]}$')
plt.ylabel(r'$\mathrm{Waist \ Location \ [cm]}$')
plt.grid(); plt.legend(); plt.show()

plt.title("Waist Values:  Propagation vs Theory")
plt.semilogy(thick_arr, betas_arr*100, label = "Waist values")
plt.semilogy(thick_arr, betap_arr*100, label = "Thin prediction")
plt.xlabel(r'$\mathrm{Lens \ Thickness \ [\mu m]}$')
plt.ylabel(r'$\mathrm{Waist \ \beta^* \ [cm]}$')
plt.grid(); plt.legend(); plt.show()

#Now for chromaticity
perbeta_arr = np.zeros(num)
perdbeta_arr = np.zeros(num)
e_spec = np.array([0, -0.01, 0.01]) + 1.0
for j in range(len(thick_arr)):
    zmult=1
    gammab = PProp.def_gamma
    tpl_n = 1.
    tpl_l = thick_arr[j]
    tpl_f = Foc.Calc_Focus_Square_CM_UM(tpl_n*1e17, tpl_l, gammab)/100
    
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
    
    argon_params = PProp.ReturnDefaultPlasmaParams(path, plasma_start = z_offset,
                                                       nset = tpl_n)
    argon_params['Z'] = z_arr[-1]*1e6
    argon_params['Nz']= len(z_arr)
    argon_params['l_flattop'] = np.NaN; argon_params['sigma_in'] = np.NaN; argon_params['sigma_out'] = np.NaN
    argon = PProp.NoPlasma_ThinPlasmaLens(argon_params, n_arr, tpl_offset*1e6 + position_error, tpl_n, tpl_l, debug)
    
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
    
    perbeta_arr[j] = 100*(betaf[1]-betaf[2])/betaf[0]
    perdbeta_arr[j] =100*(dbetaf[1]-dbetaf[2])/dbetaf[0]
    
plt.title("Chromaticity vs lens thickness at " + r'$\beta(z)=11 \ cm$')
plt.plot(thick_arr, perbeta_arr, label = r'$\Delta\beta$')
plt.plot(thick_arr, perdbeta_arr, label = r'$\Delta \partial \beta$')
plt.xlabel(r'$\mathrm{Lens \ Thickness \ [\mu m]}$')
plt.ylabel(r'$\mathrm{Percent \ Range \ [\%]}$')
plt.grid(); plt.legend(); plt.show()

plt.title("Chromaticity vs focal length at " + r'$\beta(z)=11 \ cm$')
plt.plot(focus_arr*100, perbeta_arr, label = r'$\Delta\beta$')
plt.plot(focus_arr*100, perdbeta_arr, label = r'$\Delta \partial \beta$')
plt.xlabel(r'$\mathrm{Focal \ Length \ [cm]}$')
plt.ylabel(r'$\mathrm{Percent \ Range \ [\%]}$')
plt.grid(); plt.legend(); plt.show()