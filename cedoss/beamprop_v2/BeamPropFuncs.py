#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:24:57 2018

Functions for beam propagation using Robert's code

@author: chris
"""

import sys
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.integrate as Int

sys.path.insert(0, "../../python")

import numpy as np
from beam.beams import electronbeam
from beam.elements import plasma_1d as plasma
from beam import interactions
import beam.calc.electron as ecalc
from ionization import ionization
from lens import profile
plt.rcParams['animation.ffmpeg_path'] = '/home/chris/anaconda3/envs/CU-PWFA/bin/ffmpeg'
import matplotlib.animation as animation
import scipy.constants as const
me = const.physical_constants['electron mass energy equivalent in MeV'][0]

sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim

def_startloc = 0.80
def_lenflat = 0.50 #0.50
def_nset = 0.3 #1203.7 for gas cell
def_betastar = 0.10
def_betaoffs = 0#-0.387
def_gamma = 19569.5 #10 GeV beam
def_emit = 3e-6 #New value as of July 2018
def_deltaE = 0.05#0.01#0.0025
def_sigma_hw = 14.0e4

def ReturnDefaultElectronParams(path, beta_star=def_betastar, beta_offset=def_betaoffs,
                                plasma_start=def_startloc, gamma=def_gamma, emit=def_emit,
                                deltaE = def_deltaE):
    beta_init = beta_star + np.square(plasma_start + beta_offset)/beta_star
    alpha_init = (plasma_start + beta_offset)/beta_star
    
    electronParams = {
        'name' : 'TestBeam',
        'path' : path,
        'load' : False,
        'N' : 100000,#10000,         #10000 for normal, 1000000 for production
        'gamma' : gamma,
        'emittance' : emit,
        'betax' : beta_init,
        'betay' : beta_init,
        'alphax' : alpha_init,
        'alphay' : alpha_init,
        'sigmaz' : 5e-6,
        'dE' : deltaE
    }
    return electronParams

def GaussianBeam(electronParams, debug = 0):
    beam = electronbeam.GaussianElectronBeam(electronParams)
    if debug == 1: beam.plot_current_phase();
    return beam

def VorpalBeam(path, filename, threshold, minz=-np.inf, maxz = np.inf, debug = 0):
    vorpalParams = {
            'N' : 1,
            'load' : False,
            'name' : 'VorpalBeam',
            'path' : path,
            'filename' : filename,
            'thresh' : threshold,
            'minz' : minz,
            'maxz' : maxz
    }
    beam = electronbeam.VorpalElectronBeam(vorpalParams)
    if debug ==1: beam.plot_current_phase();
    return beam

def dgammadz(ne): #Used to have old small n term
    npl0 = def_nset; npl = ne
    if npl < 1e-18:
        return 0
    else:
        dgds0 = 16.66e9 * np.sqrt(npl0/0.5) / 511e3
        #if (npl > 1/4*npl0):
        dgds = dgds0*np.sqrt(npl/npl0)*(2*np.sqrt(npl/npl0)-1)
        #else:
        #    dgds = (-dgds0)*np.sqrt(npl/npl0)*np.sin(2*np.sqrt(npl/npl0)+np.pi/2-1)
        return dgds

def dgammadz_preRobert(ne):
    npl0 = 0.5; npl = ne
    dgds0 = np.sqrt(0.5) * 1.96e-2
    if (npl > 4/9*npl0):
        dgds = dgds0*np.sqrt(npl/npl0)*(3*np.sqrt(npl/npl0)-2)
    else:
        dgds = (-dgds0)*(2/np.pi)*np.sqrt(npl/npl0)*np.sin((3*np.pi/2)*np.sqrt(npl/npl0))
    return dgds

def dgammadz_wrong(ne):
    return np.sqrt(ne) * 1.96e-2

def dgammadz_none(ne):
    return 0.0

def ReturnDefaultPlasmaParams(path, plasma_start = def_startloc,
                              nset = def_nset, sigma_hw = def_sigma_hw, scaledown = 1):
    Nx = 1;  Ny = 1
    Z = (2*plasma_start + def_lenflat)*1e6
    Nz = int((Z/scaledown)+1)
    sigma = sigma_hw/(np.sqrt(2*np.log(2)))
    plasmaParams ={
        'name' : 'TestPlasma',
        'path' : path,
        'load' : False,
        'Nx' : Nx,
        'Ny' : Ny,
        'Nz' : Nz,
        'X' : 3,
        'Y' : 3,
        'Z' : Z,
        'n0' : nset,
        'z0' : plasma_start * 1e6,
        'l_flattop' : def_lenflat*1e6,
        'sigma_in' : sigma,
        'sigma_out' : sigma,
        'atom' : ionization.Ar,
        'cyl' : False,
        'dgammadz' : dgammadz
    }
    return plasmaParams

def FineSpacingLens(z_orig, tpl_center, tpl_length):
    z_tpl = np.linspace(tpl_center-.5*tpl_length, tpl_center-.5*tpl_length, int(tpl_length))
    z = np.append(z_orig,z_tpl)
    z = np.array(sorted(z))
    return z

############# GAUSSIAN RAMPS ################################

def GaussianRampPlasma(plasmaParams, debug = 0):
    Nz = plasmaParams['Nz']
    
    nez = .5/.4995*plasmaParams['n0']*profile.plasma_gaussian_ramps(plasmaParams['z0'],
       plasmaParams['l_flattop'], plasmaParams['sigma_in'], plasmaParams['sigma_out'], Nz, plasmaParams['Z'])[1]
    
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones(Nz, dtype='double')
    argon.initialize_plasma(n, nez)
    
    if debug == 1: argon.plot_long_density_center();
    """
    Nx = plasmaParams['Nx']; Ny = plasmaParams['Ny']; Nz = plasmaParams['Nz']
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones((Nx, Ny, Nz), dtype='double')
    ne = np.zeros((Nx, Ny, Nz), dtype='double')
    ne[int(Nx/2), int(Ny/2), :] = plasmaParams['n0']*profile.plasma_gaussian_ramps(plasmaParams['z0'],
       plasmaParams['l_flattop'], plasmaParams['sigma_in'], plasmaParams['sigma_out'], Nz, plasmaParams['Z'])[1]
    
    argon.initialize_plasma(n, ne)
    if debug == 1: argon.plot_long_density_center();"""
    return argon

def GaussianRampPlasma_ThinPlasmaLens(plasmaParams, tpl_offset, tpl_n, tpl_l, debug = 0):
    tpl_l = int(tpl_l)
    Nz = plasmaParams['Nz']
    
    nez = plasmaParams['n0']*profile.plasma_gaussian_ramps(plasmaParams['z0'],
       plasmaParams['l_flattop'], plasmaParams['sigma_in'], plasmaParams['sigma_out'], Nz, plasmaParams['Z'])[1]
    
    tpl_z = plasmaParams['z0'] + tpl_offset
    dz = plasmaParams['Z']/(plasmaParams['Nz']-1)
    
    tpl_left = tpl_z-.5*tpl_l
    tpl_right = tpl_z+.5*tpl_l
    tpl_1 = int(np.floor(tpl_left/dz))
    tpl_2 = int(np.ceil(tpl_right/dz))

    lens_loc = np.array(range(tpl_2 - tpl_1 - 1)) + tpl_1 + 1
    
    for i in lens_loc:
        nez[i] = tpl_n
    
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones(Nz, dtype='double')
    argon.initialize_plasma(n, nez)
    
    if debug == 1: argon.plot_long_density_center();
    return argon

#### FOR A STUDY ON EMPIRICAL MATCHING, DIFFERENT RAMPS ####

def GaussianExitRamp(plasmaParams, sig_hw, debug = 0):
    Nz = plasmaParams['Nz']
    z_arr = np.linspace(0, plasmaParams['Z'], Nz)
    nez = np.ones(Nz, dtype='double')*plasmaParams['n0']
    sig = sig_hw/np.sqrt(2*np.log(2))
    nez = nez*np.exp(-1*np.square(z_arr)/(2*np.square(sig)))
    
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones(Nz, dtype='double')
    argon.initialize_plasma(n, nez)
    
    if debug == 1: argon.plot_long_density_center();
    return argon

################## CUSTOM OPTIONS ##########################

def CustomPlasma(plasmaParams, nez, debug = 0):
    Nz = plasmaParams['Nz']
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones(Nz, dtype='double')
    argon.initialize_plasma(n, nez)
    
    if debug == 1: argon.plot_long_density_center();
    return argon

def CustomPlasma_ThinPlasmaLens(plasmaParams, nez, tpl_offset, tpl_n, tpl_l, debug = 0):
    #tpl_l = int(tpl_l)
    Nz = plasmaParams['Nz']
    dz = plasmaParams['Z']/(plasmaParams['Nz']-1)
    tpl_z = (plasmaParams['z0'] + tpl_offset)

    tpl_left = tpl_z-.5*tpl_l
    tpl_right = tpl_z+.5*tpl_l
    tpl_1 = int(np.floor(tpl_left/dz))
    tpl_2 = int(np.ceil(tpl_right/dz))
    if debug==1: print(tpl_1,tpl_2);
    lens_loc = np.array(range(tpl_2 - tpl_1 - 1)) + tpl_1 + 1
    
    for i in lens_loc:
        nez[i] = tpl_n
    
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones(Nz, dtype='double')
    argon.initialize_plasma(n, nez)
    if debug == 1: argon.plot_long_density_center();
    return argon

#fit_tanh = [a,b,n0] - see ThreeDimensionAnalysis for more details
def CustomPlasma_ThinPlasmaLens_Tanh(plasmaParams, nez, tpl_offset, fit_tanh, debug = 0):
    tpl_l = 2*fit_tanh[0] + 6*fit_tanh[1]
    Nz = plasmaParams['Nz']
    
    tpl_z = plasmaParams['z0'] + tpl_offset
    dz = plasmaParams['Z']/(plasmaParams['Nz']-1)
    
    tpl_left = tpl_z-.5*tpl_l
    tpl_right = tpl_z+.5*tpl_l
    tpl_1 = int(np.floor(tpl_left/dz))
    tpl_2 = int(np.ceil(tpl_right/dz))
    if debug==1: print(tpl_1,tpl_2);
    lens_loc = np.array(range(tpl_2 - tpl_1 - 1)) + tpl_1 + 1
    
    z_tpl = np.linspace(-.5 * tpl_l, .5 * tpl_l, len(lens_loc))
    n_tpl = ThrDim.DoubleTanh(fit_tanh, z_tpl)
    j=0
    for i in lens_loc:
        if nez[i] < n_tpl[j]:
            nez[i] = n_tpl[j]
        j+=1
    
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones(Nz, dtype='double')
    argon.initialize_plasma(n, nez)
    if debug == 1: argon.plot_long_density_center();
    return argon

#fit_tanh = [a,b,n0] - see ThreeDimensionAnalysis for more details
#fit_lorentz = [A,gam,x0] - see ThreeDimensionAnalysis for more details
def CustomPlasma_ThinPlasmaLens_TanhLorentz(plasmaParams, nez, tpl_offset, fit_tanh, fit_lorentz, debug = 0):
    tpl_l = 2*fit_tanh[0] + 6*fit_tanh[1]
    Nz = plasmaParams['Nz']
    
    tpl_z = plasmaParams['z0'] + tpl_offset
    dz = plasmaParams['Z']/(plasmaParams['Nz']-1)
    
    tpl_left = tpl_z-.5*tpl_l
    tpl_right = tpl_z+.5*tpl_l
    tpl_1 = int(np.floor(tpl_left/dz))
    tpl_2 = int(np.ceil(tpl_right/dz))
    if debug==1: print(tpl_1,tpl_2);
    lens_loc = np.array(range(tpl_2 - tpl_1 - 1)) + tpl_1 + 1
    
    z_tpl = np.linspace(-.5 * tpl_l, .5 * tpl_l, len(lens_loc))
    n_tpl = ThrDim.DoubleTanh_Lorentzian(fit_tanh, fit_lorentz, z_tpl)
    j=0
    for i in lens_loc:
        nez[i] = n_tpl[j]
        j+=1
    
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones(Nz, dtype='double')
    argon.initialize_plasma(n, nez)
    if debug == 1: argon.plot_long_density_center();
    return argon

######################## LONE TPL #####################################

def NoPlasma_ThinPlasmaLens(plasmaParams, nez, tpl_offset, tpl_n, tpl_l, debug = 0):
    #tpl_l = int(tpl_l)
    Nz = plasmaParams['Nz']
    
    tpl_z = plasmaParams['z0'] + tpl_offset
    dz = plasmaParams['Z']/(plasmaParams['Nz']-1)
    
    tpl_left = tpl_z-.5*tpl_l
    tpl_right = tpl_z+.5*tpl_l
    tpl_1 = int(np.floor(tpl_left/dz))
    tpl_2 = int(np.ceil(tpl_right/dz))
    if debug==1: print(tpl_1,tpl_2);
    lens_loc = np.array(range(tpl_2 - tpl_1 - 1)) + tpl_1 + 1
    for i in lens_loc:
        nez[i] = tpl_n
    
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones(Nz, dtype='double')
    argon.initialize_plasma(n, nez)
    if debug == 1: argon.plot_long_density_center();
    return argon

#tpl_fit = [a,b,n0] - see ThreeDimensionAnalysis for more details
def NoPlasma_ThinPlasmaLens_Tanh(plasmaParams, nez, tpl_offset, tpl_fit, debug = 0):
    tpl_l = 2*tpl_fit[0] + 6*tpl_fit[1]
    Nz = plasmaParams['Nz']
    
    tpl_z = plasmaParams['z0'] + tpl_offset
    dz = plasmaParams['Z']/(plasmaParams['Nz']-1)
    
    tpl_left = tpl_z-.5*tpl_l
    tpl_right = tpl_z+.5*tpl_l
    tpl_1 = int(np.floor(tpl_left/dz))
    tpl_2 = int(np.ceil(tpl_right/dz))
    if debug==1: print(tpl_1,tpl_2);
    lens_loc = np.array(range(tpl_2 - tpl_1 - 1)) + tpl_1 + 1
    
    z_tpl = np.linspace(-.5 * tpl_l, .5 * tpl_l, len(lens_loc))
    n_tpl = ThrDim.DoubleTanh(tpl_fit, z_tpl)
    j=0
    for i in lens_loc:
        nez[i] = n_tpl[j]
        j+=1
    
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones(Nz, dtype='double')
    argon.initialize_plasma(n, nez)
    if debug == 1: argon.plot_long_density_center();
    return argon

########################  #####################################

def Calc_CSParams(beamParams, n_arr, z_arr):
    beta0 = beamParams['betax']
    alpha0 = beamParams['alphax']
    gb0 = beamParams['gamma']
    ne0 = max(n_arr)
    dgdz0 = 16.66e9 * np.sqrt(ne0/0.5) / 511e3
    
    return ecalc.cs_propagation(z_arr, n_arr, beta0, alpha0, gb0, dgdz0, ne0)
    #beta, alpha, gamma, gb = ecalc.cs_propagation(z, ne, beta0, alpha0, gb0, dgdz0, ne0)

def Calc_Bmag(beamParams, n_arr, z_arr):
    ne0 = max(n_arr)
    beta, alpha, gamma, gb = Calc_CSParams(beamParams, n_arr, z_arr)
    #plt.plot(z_arr,beta)
    #plt.ylim(0,0.2)
    #plt.show()
    kp = 5.95074e4 * np.sqrt(ne0)
    kb = kp/np.sqrt(2*gb[-1])
    Bmag = 0.5*(beta[-1]*kb+gamma[-1]/kb)
    return Bmag

def Calc_Proj_CSParams(beamParams, n_arr, z_arr, delta):
    gb0 = beamParams['gamma']
    delta_arr = np.linspace(-delta, delta, 101)
    gb_arr = gb0*(delta_arr+1)
    beta_arr = np.zeros(len(gb_arr))
    alpha_arr = np.zeros(len(gb_arr))
    gamma_arr = np.zeros(len(gb_arr))
    
    beta = np.zeros((len(gb_arr),len(z_arr)))
    alpha = np.zeros(beta.shape)
    gamma = np.zeros(beta.shape)
    
    for n in range(len(gb_arr)):
        beamParamsCopy = beamParams.copy()
        beamParamsCopy['gamma'] = gb_arr[n]
        beta[n], alpha[n], gamma[n], gb = Calc_CSParams(beamParamsCopy, n_arr, z_arr)
    
    bmag_arr = np.zeros(len(z_arr))
    beta_pro_arr = np.zeros(len(z_arr))
    for j in range(len(bmag_arr)):
        for i in range(len(gb_arr)):
            index = j 
            beta_arr[i] = beta[i][index]; alpha_arr[i] = alpha[i][index]; gamma_arr[i] = gamma[i][index]

        #beta_pro = Int.simps(beta_arr, delta_arr)/(delta_arr[-1]-delta_arr[0])
        #alpha_pro = Int.simps(alpha_arr, delta_arr)/(delta_arr[-1]-delta_arr[0])
        #gamma_pro = Int.simps(gamma_arr, delta_arr)/(delta_arr[-1]-delta_arr[0])
        beta_pro = np.average(beta_arr)
        alpha_pro = np.average(alpha_arr)
        gamma_pro = np.average(gamma_arr)
        
        beta_pro_arr[j] = beta_pro
        
        bmag_arr[j] = np.sqrt(beta_pro*gamma_pro - alpha_pro**2)
        
    return gb_arr, beta_arr, alpha_arr, gamma_arr, bmag_arr, beta_pro_arr

def Plot_CSEvo(beamParams, n_arr, z_arr, z_offset = 0, legend_loc=0, subset = False):
    beta, alpha, gamma, gb = Calc_CSParams(beamParams, n_arr, z_arr)
    beta0, alpha0, gamma0, gb0 = Calc_CSParams(beamParams, np.zeros(len(z_arr)), z_arr)
    
    z_arr = z_arr - z_offset
    #print(" ","beta","beta0"); print(" ",beta[740],beta0[740])
    
    if subset is not False:
        #subset = [start, end]
        z_arr = z_arr[subset[0]:subset[1]]
        n_arr = n_arr[subset[0]:subset[1]]
        beta = beta[subset[0]:subset[1]]
        beta0 = beta0[subset[0]:subset[1]]
        
    fig, ax1 = plt.subplots(figsize=(6,4))
    #plt.title("Beta function evolution at "+r'$n_0=$'+str(max(n_arr))+r'$\,\mathrm{\times 10^{17}cm^{-3}}$')
    ax1.plot(z_arr*1e2, np.array(beta)*1e2, 'b-', label=r'$\beta$')
    ax1.plot(z_arr*1e2, np.array(beta0)*1e2, 'b--',label=r'$\beta_{vac}$')
    ax1.set_ylabel(r'$\beta\,\mathrm{(cm)}$', color = 'b')
    ax1.tick_params('y', colors = 'b')
    ax1.set_xlabel('z (cm)')
    ax1.set_ylim([-0.05,20.05])
    #ax1.set_ylim([4.8,5.4])
    
    ax2 = ax1.twinx()
    ax2.plot(z_arr*1e2, n_arr/max(n_arr), 'g-')
    ax2.plot(z_arr*1e2, n_arr/max(n_arr), 'g-')
    ax2.set_ylabel(r'$n/n_0$',color = 'g')
    ax2.tick_params('y', colors = 'g')
    ###
    #ax1.set_xlim([-8.5,-6.9])
    #ax2.set_xlim([-8.5,-6.9])
    #ax2.set_ylim([1e-4,2])
    ###
    
    ax1.grid()
    ax1.legend(loc=legend_loc); plt.tight_layout()
    #plt.savefig('/home/chris/Desktop/fig3.eps',format='eps',bbox_inches='tight',dpi=150)
    plt.show()
    
    return beta,alpha,gamma,gb
    
#Changed this a little bit, Can get rid of the 'p' lines 
def Plot_CSEvo_MatchedCompare(beamParams, beamParams_matched, n_arr, z_arr, z_offset = 0, legend_loc=0):
    n_arrp = np.copy(n_arr)
    n_arrp[:int(len(n_arrp)/3)] = 0
    
    beta, alpha, gamma, gb = Calc_CSParams(beamParams, n_arr, z_arr)
    beta0, alpha0, gamma0, gb0 = Calc_CSParams(beamParams, np.zeros(len(z_arr)), z_arr)
    beta0p, alpha0p, gamma0p, gb0p = Calc_CSParams(beamParams, n_arrp, z_arr)
    betaM, alphaM, gammaM, gbM = Calc_CSParams(beamParams_matched, np.zeros(len(z_arr)), z_arr)
    
    z_arr = z_arr - z_offset
    #print(" ","beta","beta0"); print(" ",beta[740],beta0[740])
    fig, ax1 = plt.subplots(figsize=(13,5))
    plt.rcParams.update({'font.size': 12})
    #plt.title("Beta function evolution at "+r'$n_0=$'+str(max(n_arr))+r'$\,\mathrm{\times 10^{17}cm^{-3}}$')
    ax1.semilogy(z_arr*1e2, np.array(beta0p)*1e2, 'r-',label=r'$\beta_{mismatch}$')
    ax1.plot(z_arr*1e2, np.array(beta)*1e2, 'b-', label=r'$\beta$')
    ax1.plot(z_arr*1e2, np.array(beta0)*1e2, 'r--',label=r'$\beta_{vac}$')
    ax1.plot(z_arr*1e2, np.array(betaM)*1e2, 'b--',label=r'$\beta_{m,vac}$')
    ax1.set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'b')
    ax1.tick_params('y', colors = 'b')
    ax1.set_xlabel('z [cm]')
    ax1.set_ylim([-0.05,1400.05])
    
    ax2 = ax1.twinx()
    ax2.plot(z_arr*1e2, n_arr/max(n_arr), 'g-',lw=2.5)
    ax2.set_ylabel(r'$n/n_0$',color = 'g')
    ax2.tick_params('y', colors = 'g')
    ax1.grid(); ax1.legend(loc=legend_loc); plt.show()

def Plot_CSEvo_FinalCompare(beamParams, n_arr, z_arr, z_offset = 0, legend_loc=0, subset = False, plot = 1):
    beta, alpha, gamma, gb = Calc_CSParams(beamParams, n_arr, z_arr)
    
    beamParams0 = beamParams.copy()
    betastar_fin = 1/gamma[-1]
    waist_fin = alpha[-1]*betastar_fin + z_arr[-1]
    alpha_fin = waist_fin/betastar_fin
    beta_fin = betastar_fin + waist_fin**2/betastar_fin
    beamParams0['alphax'] = alpha_fin;  beamParams0['alphay'] = alpha_fin
    beamParams0['betax'] = beta_fin;  beamParams0['betay'] = beta_fin
    beta0, alpha0, gamma0, gb0 = Calc_CSParams(beamParams0, np.zeros(len(z_arr)), z_arr)
    
    z_arr = z_arr - z_offset
    #print(" ","beta","beta0"); print(" ",beta[740],beta0[740])
    
    if subset is not False:
        #subset = [start, end]
        z_arr = z_arr[subset[0]:subset[1]]
        n_arr = n_arr[subset[0]:subset[1]]
        beta = beta[subset[0]:subset[1]]
        beta0 = beta0[subset[0]:subset[1]]
    
    if plot == 1:
        fig, ax1 = plt.subplots()
        plt.title("Beta function evolution at "+r'$n_0=$'+str(max(n_arr))+r'$\,\mathrm{\times 10^{17}cm^{-3}}$')
        ax1.plot(z_arr*1e2, np.array(beta)*1e2, 'b-', label=r'$\beta$')
        ax1.plot(z_arr*1e2, np.array(beta0)*1e2, 'b--',label=r'$\beta_{vac}$')
        ax1.set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'b')
        ax1.tick_params('y', colors = 'b')
        ax1.set_xlabel('z [cm]')
        ax1.set_ylim([-0.05,20.05])
        
        ax2 = ax1.twinx()
        ax2.plot(z_arr*1e2, n_arr/max(n_arr), 'g-')
        ax2.set_ylabel(r'$n/n_0$',color = 'g')
        ax2.tick_params('y', colors = 'g')
        ax1.grid(); ax1.legend(loc=legend_loc); plt.show()
    
    return beta,alpha,gamma,gb

def PropBeamPlasma(beam, plasma, z_arr, dumpPeriod, cores=4, debug = 0):
    m = int(len(z_arr)/dumpPeriod) -1
    interactions.electron_plasma(beam, plasma, z_arr, dumpPeriod, cores)
    if debug == 1:
        beam.plot_current_phase()
        print('Initial emittance:', np.average(beam.get_emittance_n(0))*1e6, 'mm.mrad')
        print('Final emittance:', np.average(beam.get_emittance_n(m))*1e6, 'mm.mrad')
    return beam

def GetBmag(beam,m):
    e0 = np.average(beam.get_emittance_n(0))*1e6
    em = np.average(beam.get_emittance_n(m))*1e6
    return em/e0

def PlotEmittance(beam, z_arr, m):
    en = np.zeros(m, dtype='double')
    s = np.zeros(m, dtype='double')
    for i in range(m):
        #j = int(i * len(z_arr)/m)
        en[i] = np.average(beam.get_emittance_n(i))*1e6
        s[i] = beam.get_save_z(i)
        #s[i] = z_arr[j]*1e-4
    plt.plot(s, en)
    plt.title("Emittance Evolution")
    plt.xlabel("s [cm]")
    plt.ylabel("normalized emittance [mm-mrad]")
    plt.grid(); plt.show()
    print(en[0])
    return

def GetSigmaMin(beam, m):
    sig = np.zeros(m, dtype='double')
    for i in range(m):
        sig[i] = np.average(beam.get_sigmar(i))
    return min(sig)

def GetBetaMin(beam, m):
    sig = np.zeros(m, dtype='double')
    for i in range(m):
        sig[i] = np.average(beam.get_sigmar(i))
    ef = np.average(beam.get_emittance(m))
    return np.square(min(sig))/ef

def PlotSigmar(beam, z_arr, m):
    sig = np.zeros(m, dtype='double')
    sigx = np.zeros(m, dtype='double')
    sigy = np.zeros(m, dtype='double')
    s = np.zeros(m, dtype='double')
    for i in range(m):
        #j = int(i * len(z_arr)/m)
        sig[i] = np.average(beam.get_sigmar(i))*1e6
        sigx[i] = beam.get_sigmar(i)[0]*1e6
        sigy[i] = beam.get_sigmar(i)[1]*1e6
        s[i] = beam.get_save_z(i)*1e2
        #s[i] = z_arr[j]*1e-4
        
    plt.plot(s, sig, label='average')
    plt.plot(s, sigx, label='x')
    plt.plot(s, sigy, label='y')
    #plt.title("Sigmar Evolution")
    plt.xlabel(r'$z \ \mathrm{(cm)}$')
    plt.ylabel(r'$\sigma_{x,y} \ \mathrm{(\mu m)}$')
    plt.grid(); plt.legend(); plt.show()
    return sig, s

def PlotGamma(beam, z_arr, m):
    gam = np.zeros(m, dtype='double')
    s = np.zeros(m, dtype='double')
    for i in range(m):
        j = int(i * len(z_arr)/m)
        gam[i] = np.average(beam.get_gamma_n(i))
        s[i] = z_arr[j]*1e-4
    
    plt.plot(s, gam)
    plt.title("Gamma Evolution")
    plt.xlabel("s [cm]")
    plt.ylabel("gamma")
    plt.grid(); plt.show()
    print(gam[-1])
    return

def GetFinalGamma(beam, m):
    return np.average(beam.get_gamma_n(m))

def PlotEmittance_Compare(beam, beam2, z_arr, m, first, second):
    en = np.zeros(m, dtype='double')
    en2 = np.zeros(m, dtype='double')
    s = np.zeros(m, dtype='double')
    for i in range(m):
        j = int(i * len(z_arr)/m)
        en[i] = np.average(beam.get_emittance_n(i))*1e6
        en2[i] = np.average(beam2.get_emittance_n(i))*1e6
        s[i] = z_arr[j]*1e-4
    
    plt.plot(s, en,label = first)
    plt.plot(s, en2, label = second)
    plt.title("Emittance Evolution")
    plt.xlabel("s [cm]")
    plt.ylabel("normalized emittance [mm-mrad]")
    plt.grid(); plt.legend(); plt.show()
    return

def PlotSigmar_Compare(beam, beam2, z_arr, m, first, second):
    sig = np.zeros(m, dtype='double')
    sig2 = np.zeros(m, dtype='double')
    s = np.zeros(m, dtype='double')
    for i in range(m):
        j = int(i * len(z_arr)/m)
        sig[i] = np.average(beam.get_sigmar(i))*1e6
        sig2[i] = np.average(beam2.get_sigmar(i))*1e6
        s[i] = z_arr[j]*1e-4
    
    plt.plot(s, sig, label = first)
    plt.plot(s, sig2, label = second)
    plt.title("Sigmar Evolution")
    plt.xlabel("s [cm]")
    plt.ylabel("sigmar [um]")
    plt.grid(); plt.legend(); plt.show()
    return

def PlotContour(contour, x_arr, y_arr, x_label, y_label, simple = False, swapx = 0, swapy = 0):
        # find location of min(B)
        
    y_arr = y_arr-swapy
    x_arr = x_arr-swapx
    
    i_Bmin_x = np.argmin(np.min(contour,1))
    i_Bmin_y = np.argmin(np.min(contour,0))
    Bmin_x = x_arr[i_Bmin_x]
    Bmin_y = y_arr[i_Bmin_y]
    Bmin   = np.min(np.min(contour))
    
    print('matching x value: ',Bmin_x)
    print('matching y value: ',Bmin_y)
    print('matched B-mag: ',Bmin)
    print('emittance growth (%): ',100*(Bmin-1)/Bmin)
    
    X = np.tile(x_arr.reshape(-1,1),(1,len(y_arr)))
    Y = np.tile(y_arr.T,(len(x_arr),1))
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.contourf(X,Y,np.log10(contour),100,\
                cmap=plt.get_cmap('Vega20c'),\
                linewidth=2.0)
    cbar = plt.colorbar()
    plt.scatter(Bmin_x,Bmin_y,color='k')
    cbar.ax.set_ylabel(r'$log_{10}(B_m)$')
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.show()
#    plt.title(r'beam matching for %s ramp'%shape_up)

    # thin line contour map of B
    
    plt.figure(figsize=(7,5))
    plt.rcParams.update({'font.size': 12})
    
    levels = np.array([1.01,1.05,1.1,1.2,1.5,2.0,3.0,4.0,5.0])
    labels = np.array([1.01,1.05,1.1,1.2,1.5,2.0,3.0,4.0,5.0])
    fig, axes = plt.subplots(1,1, sharey=True)
    CS = plt.contour(X,Y,contour,levels,cmap=plt.get_cmap('Vega20b'))
    plt.clabel(CS,labels,fontsize=10,inline=1,fmt='%1.2f')
    if simple:
        plt.grid()
        print()
    else:
        cbar = plt.colorbar()
    plt.scatter(Bmin_x-swapx,Bmin_y-swapy,color='k')
    cbar.ax.set_ylabel(r'$B_m$')
    cbar.set_ticks(levels)
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.tight_layout()
    #plt.savefig('/home/chris/Desktop/fig2.eps',format='eps',bbox_inches='tight',dpi=150)
    plt.show()
#    plt.title(r'beam matching for %s ramp'%shape_up)

    return [Bmin_x,Bmin_y]

def PlotContour_General(contour, x_arr, y_arr, x_label, y_label, data_label):
        # find location of min(B)
    i_Bmin_x = np.argmin(np.min(contour,1))
    i_Bmin_y = np.argmin(np.min(contour,0))
    Bmin_x = x_arr[i_Bmin_x]
    Bmin_y = y_arr[i_Bmin_y]
    Bmin   = np.min(np.min(contour))
    
    print('matching x value: ',Bmin_x)
    print('matching y value: ',Bmin_y)
    print('minimum value: ',Bmin)
    
    X = np.tile(x_arr.reshape(-1,1),(1,len(y_arr)))
    Y = np.tile(y_arr.T,(len(x_arr),1))
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.contourf(X,Y,contour,100,\
                cmap=plt.get_cmap('Vega20c'),\
                linewidth=2.0)
    cbar = plt.colorbar()
    plt.scatter(Bmin_x,Bmin_y,color='k')
    cbar.ax.set_ylabel(data_label)
    plt.ylabel(y_label)
    plt.xlabel(x_label)

    return [Bmin_x,Bmin_y]

def Plot_CSEvo_MatchedCompare_NicePlot(beamParams, beamParams_matched, n_arr, z_arr, z_offset = 0, legend_loc=0):
    beta, alpha, gamma, gb = Calc_CSParams(beamParams, n_arr, z_arr)
    beta0, alpha0, gamma0, gb0 = Calc_CSParams(beamParams, np.zeros(len(z_arr)), z_arr)
    betaM, alphaM, gammaM, gbM = Calc_CSParams(beamParams_matched, np.zeros(len(z_arr)), z_arr)
    
    z_arr = z_arr - z_offset
    #print(" ","beta","beta0"); print(" ",beta[740],beta0[740])
    """
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}

    plt.rc('font', **font)
    """
    lwid = 2.0
    
    fig, ax1 = plt.subplots()
    ax1.plot(z_arr*1e2, np.array(beta)*1e2, 'b-', label=r'$\beta$', linewidth = lwid)
    ax1.plot(z_arr*1e2, np.array(beta0)*1e2, 'b--',label=r'$\beta_{vac}$', linewidth = lwid)
    ax1.plot(z_arr*1e2, np.array(betaM)*1e2, 'r--',label=r'$\beta_{m,vac}$', linewidth = lwid)
    ax1.set_ylabel(r'$\beta\,\mathrm{[cm]}$', color = 'b')
    ax1.tick_params('y', colors = 'b')
    ax1.set_xlabel('z [cm]')
    
    ax2 = ax1.twinx()
    ax2.plot(z_arr*1e2, n_arr/max(n_arr), 'g-', linewidth = lwid)
    ax2.set_ylabel(r'$n/n_0$',color = 'g')
    ax2.tick_params('y', colors = 'g')
    ax1.set_xlim([-45, 5])
    ax1.set_ylim([0, 20])
    ax2.set_ylim([0, 1.1])
    fig.set_size_inches(3.4*1.5,2.5*1.5)
    ax1.grid(); ax1.legend(); plt.show()

###############################################################################

def GaussianRampPlasma_ThinPlasmaLens_Disused(plasmaParams, tpl_offset, tpl_n, tpl_l):
    tpl_l = int(tpl_l)
    Nx = plasmaParams['Nx']; Ny = plasmaParams['Ny']; Nz = plasmaParams['Nz']
    
    nez = plasmaParams['n0']*profile.plasma_gaussian_ramps(plasmaParams['z0'],
       plasmaParams['l_flattop'], plasmaParams['sigma_in'], plasmaParams['sigma_out'], Nz, plasmaParams['Z'])[1]
    
    tpl_z = plasmaParams['z0'] + tpl_offset
    dz = plasmaParams['Z']/(plasmaParams['Nz']-1)
    
    tpl_left = tpl_z-.5*tpl_l
    tpl_right = tpl_z+.5*tpl_l
    tpl_1 = int(np.floor(tpl_left/dz))
    tpl_2 = int(np.ceil(tpl_right/dz))
    
    z_orig = np.linspace(0,plasmaParams['Z'],plasmaParams['Nz'])
    z_fine = FineSpacingLens(z_orig, tpl_z, tpl_l)
    
    nez_tpl = np.ones(tpl_l)*tpl_n
    
    nez_A = nez[0 : tpl_1 + 1]
    nez_B = nez[tpl_2 :]
    nez_C = nez[tpl_1 + 1 : tpl_2]
    for i in range(len(nez_C)):
        nez_C[i] = tpl_n
    
    nez_fine = np.append(np.append(nez_A, nez_tpl), np.append(nez_C, nez_B))
    plt.plot(z_fine,nez_fine)
    
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones((Nx, Ny, Nz + tpl_l), dtype='double')
    ne = np.zeros((Nx, Ny, Nz + tpl_l), dtype='double')
    ne[int(Nx/2), int(Ny/2), :] = nez_fine
    argon.initialize_plasma(n, ne)
    argon.plot_long_density_center()
    #return argon
    return argon