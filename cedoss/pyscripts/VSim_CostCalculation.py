#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 12:35:55 2019

While the VSim_Estimation scripts are useful for telling you what your parameters
should be, they often would not be totally identical to what the actual simulation
is.  For instance, the estimations would calculate dx and give a Nx, but VSim is
the opposite and takes an Nx and gives dx.

This script simply uses the direct input parameters to VSim and calculates the
expected MPP hour usage.  I will also include future parameters that scale this
"ideal expected" cost to the actual cost on NERSC to reflect times which are
harder to calculate, such as setup, dumping, and overall processor speed.

@author: chris
"""
import numpy as np

Nmal = 8
Nlong = 360;  Llong = 202*1e-6
Ntran = 518;  Ltran = 300*1e-6

SimEnd = 0.06*1e-2;  SimStart = -0.07*1e-2;  waistloc = 0
Ndumps = 10;

DriveSigR = 5.1 *1e-6;  DriveSigZ = 5.2 *1e-6;  DriveEmit = 5.3 *1e-6
WitnsSigR = 3.9 *1e-6;  WitnsSigZ = 12. *1e-6;  WitnsEmit = 3.0 *1e-6

gbC = 19569.5

## Constants

cdtfac = 0.99999
lightspeed = 299792458.0

ppc_p = 1
ppc_d = 8; ppc_w = 8

## Calculations which mirror .pre

Ntran_tot = Ntran + 2*Nmal

Dy = Ltran / Ntran_tot; Dz = Dy
Dx = Llong / Nlong

L_sim = SimEnd - SimStart + Llong

dl = 1/np.sqrt(1/(Dy)**2+1/(Dz)**2+1/(Dx)**2)
dt = cdtfac*dl/lightspeed

Lt = L_sim / lightspeed
dumpperiod = int(Lt/dt/Ndumps)
Nt = Ndumps * dumpperiod

## Calculations from Estimation scripts
#Particle Steps
Ncells = Ntran_tot**2*Nlong
Nprocs = Ncells / 29000
plasma_steps = Ncells*Nt*ppc_p

#Below is to account for the beam being larger when initialized away from waist
eff_fac_d = np.sqrt(1+np.square(DriveEmit/gbC*(SimStart-waistloc)/np.square(DriveSigR)))
eff_fac_w = np.sqrt(1+np.square(WitnsEmit/gbC*(SimStart-waistloc)/np.square(WitnsSigR)))

drive_Ntran = 2*(5*DriveSigR)/Dy*eff_fac_d
drive_Nlong = 2*(5*DriveSigZ)/Dx
drive_steps = drive_Ntran**2*drive_Nlong*Nt*ppc_d

witness_Ntran = 2*(5*WitnsSigR)/Dy*eff_fac_w
witness_Nlong = 2*(5*WitnsSigZ)/Dx
witness_steps = witness_Ntran**2*witness_Nlong*Nt*ppc_w

part_steps = plasma_steps + drive_steps + witness_steps

#MPP Core Hours
MPP = part_steps*400/3.6e12 * 2

#GB per Dump - scales w.r.t. some 3D Robert sim
size_E = (455)/1000 * (Ncells)/(250**3)
size_rho = 3*(150)/1000 * (Ncells)/(250**3)
size_drive = (2*455)/1000 * (drive_Ntran**2*drive_Nlong)/(250**3) * ppc_d
size_witness = (2*455)/1000 * (witness_Ntran**2*witness_Nlong)/(250**3) * ppc_w
size_per_dump = size_E + size_rho + size_drive + size_witness

#Total Size
GBsize = Ndumps * size_per_dump

machine = 'Haswell'
corespn = 32
nodes = np.ceil(Nprocs/corespn)

print("Transverse    ","Ny: ",str(Ntran_tot)," dy: ",str(Dy))
print("Longitudinal  ","Nx: ",str(Nlong)," dx: ",str(Dx))
print("Time          ","Nt: ",str(Nt)," dt: ",str(dt))
print("Sim Length    ","L:  ",str(L_sim), " Dumps: ",str(Ndumps))
print("# Cells, Procs","N:",str(Ncells)," Pr: ",str(Nprocs))
print("Machine       ",machine," # Nodes: ",str(nodes))
print("Storage       ","GB per dump: ",str(size_per_dump))
print("              ","Total TB: ",str(GBsize/1000))
print("MPP Hours: ",str(MPP))