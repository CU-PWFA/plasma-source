#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 12:01:26 2018

Coming up with parameter space for thin plasma lens simulation

Note this is for 2D, so this is different than previous estimation script

@author: chris
"""

import numpy as np

dens = 1e18 #cm^-3
gbC = 19569.5
ppc_p = 1
nmal = 12

z_start = -0.001
z_waist = 0.0
z_end = 0.006
L_sim = z_end - z_start #m

specie = ["drive", "witness"]
witness_delay = 50. #um
drive_window = 10. #um
lrms = [5.2, 18.] #um
sigrms = [4.1, 4.7] #um
emit_n = [5.3, 7.] #um-rad
ppc = [8,8] #and plasma is 1

#Transverse Domain
umtrans = 100 #um
k_b0 = 1.33e-4*np.sqrt(dens/gbC)
sigmar_m = [0,0]
for i in range(len(specie)):
    sigmar_m[i] = np.sqrt(emit_n[i]/1e6/gbC/k_b0)
dx = min(sigmar_m)/5
Nx = np.ceil(umtrans/1e6/dx)+2*nmal

#Longitudinal Domain

umlong = 5*lrms[0]+5*lrms[1] + witness_delay + drive_window
#dz = min(lrms)/5/1e6
dz = dx
Nz = np.ceil(umlong/1e6/dz)

#Simulation Length
dt = min([dx,dz])/3e8
Nt = np.ceil(L_sim/3e8/dt)

#Particle Steps
plasma_steps = Nx*Nz*Nt*ppc_p

eff_fac = [0,0]
for i in range(len(specie)):
    eff_fac[i] = np.sqrt(1+np.square(emit_n[i]/1e6/gbC*(z_start-z_waist)/np.square(sigrms[i]/1e6)))
print("eff_fac: ",str(eff_fac)); print()

drive_Nx = 2*(5*sigrms[0])/1e6/dx*eff_fac[0]
drive_Nz = 2*(5*lrms[0])/1e6/dz
drive_steps = drive_Nx*drive_Nz*Nt*ppc[0]

witness_Nx = 2*(5*sigrms[1])/1e6/dx*eff_fac[1]
witness_Nz = 2*(5*lrms[1])/1e6/dz
witness_steps = witness_Nx*witness_Nz*Nt*ppc[1]

part_steps = plasma_steps + drive_steps + witness_steps

#MPP Core Hours
MPP = part_steps*400/3.6e12
#NEED TO REDO FOR 2D.  DEFINITELY NOT 4 TB (OR SO I WOULD HOPE)
#GB per Dump - scales w.r.t. some 3D Robert sim
size_E = (455)/1000 * (Nx*Nz)/(250**2)
size_rho = 3*(150)/1000 * (Nx*Nz)/(250**2)
size_drive = (2*455)/1000 * (drive_Nx*drive_Nz)/(250**2) * ppc[0]
size_witness = (2*455)/1000 * (witness_Nx*witness_Nz)/(250**2) * ppc[1]
size_per_dump = size_E + size_rho + size_drive + size_witness

#Total Size
num_dumps = np.ceil(L_sim/(1/(k_b0)/5))
GBsize = num_dumps * size_per_dump

print("Transverse    ","Nx: ",str(Nx-2*nmal)," dx: ",str(dx))
print("Longitudinal  ","Nz: ",str(Nz)," dz: ",str(dz))
print("Time          ","Nt: ",str(Nt)," dt: ",str(dt))
print("Sim Length    ","L:  ",str(L_sim), " Dumps: ",str(num_dumps))
print("Storage       ","GB per dump: ",str(size_per_dump))
print("              ","Total TB: ",str(GBsize/1000))
print("MPP Hours: ",str(MPP))