#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 11:01:35 2018

For estimating a 3D TPL run

NOT YET IMPLEMENTED FOR ASYMMETRIC BEAMS

@author: chris
"""

import numpy as np

dens = 3e16 #cm^-3
gbC = 19569.5
ppc_p = 1
lightspeed = 299792458.0

z_start = -500e-6
z_waist = 0
z_end = 400e-6


specie = ["drive", "witness"]
witness_delay = 150. #um
drive_window = 26. #um
lrms = [5.2, 5.2] #um
sigrms = [10, 2.81] #um
emit_n = [3.2, 3.15] #um-rad
ppc = [8,8] #and plasma is 1

#Transverse Domain
umtrans = 300 #um
#k_b0 = 1.33e-4*np.sqrt(dens/gbC)
#sigmar_m = [0,0]
#for i in range(len(specie)):
#    sigmar_m[i] = np.sqrt(emit_n[i]/1e6/gbC/k_b0)
dx = min(sigrms)/5
Nx = np.ceil(umtrans/dx)

#Longitudinal Domain

umlong = 5*lrms[1] + witness_delay + drive_window#+5*lrms[0]
#or
#umlong = 218
L_sim = z_end - z_start + umlong*1e-6 #m

#dz = min(lrms)/5/1e6
dz = dx
Nz = np.ceil(umlong/dz)

#Simulation Length
dl = 1/np.sqrt(1/(dx*1e-6)**2+1/(dx*1e-6)**2+1/(dz*1e-6)**2)
dt = 0.99999*dl/lightspeed
Nt = np.ceil(L_sim/lightspeed/dt)

#Particle Steps
plasma_steps = Nx**2*Nz*Nt*ppc_p

eff_fac = [0,0]
for i in range(len(specie)):
    eff_fac[i] = np.sqrt(1+np.square(emit_n[i]/1e6/gbC*(z_start-z_waist)/np.square(sigrms[i]/1e6)))
print("eff_fac: ",str(eff_fac)); print()

drive_Nx = 2*(5*sigrms[0])/1e6/dx*eff_fac[0]
drive_Nz = 2*(5*lrms[0])/1e6/dz
drive_steps = drive_Nx**2*drive_Nz*Nt*ppc[0]

witness_Nx = 2*(5*sigrms[1])/1e6/dx*eff_fac[1]
witness_Nz = 2*(5*lrms[1])/1e6/dz
witness_steps = witness_Nx**2*witness_Nz*Nt*ppc[1]

part_steps = plasma_steps + drive_steps + witness_steps

#MPP Core Hours
MPP = part_steps*400/3.6e12

#GB per Dump - scales w.r.t. some 3D Robert sim
size_E = (455)/1000 * (Nx**2*Nz)/(250**3)
size_rho = 3*(150)/1000 * (Nx**2*Nz)/(250**3)
size_drive = (2*455)/1000 * (drive_Nx**2*drive_Nz)/(250**3) * ppc[0]
size_witness = (2*455)/1000 * (witness_Nx**2*witness_Nz)/(250**3) * ppc[1]
size_per_dump = size_E + size_rho + size_drive + size_witness

#Total Size
num_dumps = 10 #np.ceil(L_sim/(1/(k_b0)/5))
GBsize = num_dumps * size_per_dump

Ncells = Nx*Nx*Nz
Nprocs = Ncells/29000

print("Transverse    ","Nx: ",str(Nx)," dx: ",str(dx))
print("Longitudinal  ","Nz: ",str(Nz)," dz: ",str(dz)," Lz: ",str(umlong))
print("Time          ","Nt: ",str(Nt)," dt: ",str(dt))
print("Sim Length    ","L:  ",str(L_sim), " Dumps: ",str(num_dumps))
print("# Cells, Procs","N:",str(Ncells)," Pr:",str(Nprocs))
print("Storage       ","GB per dump: ",str(size_per_dump))
print("              ","Total TB: ",str(GBsize/1000))
print("MPP Hours: ",str(MPP))
