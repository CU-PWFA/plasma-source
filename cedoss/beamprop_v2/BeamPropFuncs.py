#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:24:57 2018

@author: chris
"""

import sys
import matplotlib.pyplot as plt
import scipy.stats as stats

sys.path.insert(0, "../../python")

import numpy as np
from beam.beams import electronbeam
from beam.elements import plasma_1d as plasma
from beam import interactions
from ionization import ionization
from lens import profile
plt.rcParams['animation.ffmpeg_path'] = '/home/chris/anaconda3/envs/CU-PWFA/bin/ffmpeg'
import matplotlib.animation as animation

plasma_start_loc = 0.75

def ReturnDefaultElectronParams(path):
    beta_star = 0.10
    beta_offset = -0.36
    plasma_start = plasma_start_loc
    
    beta_init = beta_star + np.square(plasma_start + beta_offset)/beta_star
    alpha_init = (plasma_start + beta_offset)/beta_star
    
    electronParams = {
        'name' : 'TestBeam',
        'path' : path,
        'load' : False,
        'N' : 10000,
        'gamma' : 20000,
        'emittance' : 5e-6,
        'betax' : beta_init,
        'betay' : beta_init,
        'alphax' : alpha_init,
        'alphay' : alpha_init,
        'sigmaz' : 5e-6,
        'dE' : 0.01
    }
    return electronParams

def GaussianBeam(electronParams, debug = 0):
    beam = electronbeam.GaussianElectronBeam(electronParams)
    if debug == 1: beam.plot_current_phase();
    return beam

def dgammadz(ne):
    return np.sqrt(ne) * 1.96e-2

def dgammadz_basic(ne):
    return 0.0

def ReturnDefaultPlasmaParams(path):
    Nx = 1;  Ny = 1
    Z = 2e6
    Nz = int((Z/10)+1)
    plasmaParams ={
        'name' : 'TestPlasma',
        'path' : path,
        'load' : False,
        'Nx' : Nx,
        'Ny' : Ny,
        'Nz' : Nz,
        'X' : 100,
        'Y' : 100,
        'Z' : Z,
        'n0' : 0.5,
        'z0' : plasma_start_loc * 1e6,
        'l_flattop' : 0.5e6,
        'sigma_in' : 13.25e4,
        'sigma_out' : 13.25e4,
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

def GaussianRampPlasma(plasmaParams):
    Nx = plasmaParams['Nx']; Ny = plasmaParams['Ny']; Nz = plasmaParams['Nz']
    argon = plasma.Plasma(plasmaParams)
    n = plasmaParams['n0']*np.ones((Nx, Ny, Nz), dtype='double')
    ne = np.zeros((Nx, Ny, Nz), dtype='double')
    ne[int(Nx/2), int(Ny/2), :] = plasmaParams['n0']*profile.plasma_gaussian_ramps(plasmaParams['z0'],
       plasmaParams['l_flattop'], plasmaParams['sigma_in'], plasmaParams['sigma_out'], Nz, plasmaParams['Z'])[1]
    
    argon.initialize_plasma(n, ne)
    argon.plot_long_density_center()
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

def PlotEmittance(beam,m):
    en = np.zeros(m, dtype='double')
    for i in range(m):
        en[i] = np.average(beam.get_emittance_n(i))*1e6
    
    plt.plot(en)
    plt.show()
    return

def PlotContour(contour, x_arr, y_arr, x_label, y_label):
        # find location of min(B)
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
#    plt.title(r'beam matching for %s ramp'%shape_up)

    # thin line contour map of B
    levels = np.array([1.01,1.05,1.1,1.2,1.5,2.0,3.0,4.0,5.0])
    labels = np.array([1.01,1.05,1.1,1.2,1.5,2.0,3.0,4.0,5.0])
    fig, axes = plt.subplots(1,1, sharey=True)
    CS = plt.contour(X,Y,contour,levels,cmap=plt.get_cmap('Vega20b'))
    plt.clabel(CS,labels,fontsize=9,inline=1,fmt='%1.2f')
    cbar = plt.colorbar()
    plt.scatter(Bmin_x,Bmin_y,color='k')
    cbar.ax.set_ylabel(r'$B_m$')
    cbar.set_ticks(levels)
    plt.ylabel(y_label)
    plt.xlabel(x_label)
#    plt.title(r'beam matching for %s ramp'%shape_up)

    return [Bmin_x,Bmin_y]

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