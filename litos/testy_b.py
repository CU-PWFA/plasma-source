#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 15:10:13 2017

@author: mike
"""

import numpy as np
import matplotlib.pyplot as plt
import nat_consts as nc
import particle_beam as pb
import plasma_source as ps

# define plasma bulk (flat-top) properties
npl0   = 1e17 # cm^-3
dEds0  = 20.00e9 # eV/m
dgds0  = dEds0/nc.me
L_ft   = 0.50 # m

# define plasma up-ramp
shape_up = 'gauss'
hw_up    = 0.15 # m
L_up     = 0.30 # m
top_up   = L_up # m

# define plasma down-ramp
shape_dn = shape_up
hw_dn    = hw_up # m
L_dn     = L_up # m
top_dn   = 0  # m

# define longitudinal steps
ds   = 0.001 # m
s_ft = np.linspace(0,L_ft,round(L_ft/ds+1))
s_up = np.linspace(0,L_up,round(L_up/ds+1))
s_dn = np.linspace(0,L_dn,round(L_dn/ds+1))

# make plasma
bulk    = ps.make_bulk(s_ft,npl0,dgds0)
up_ramp = ps.make_ramp(s_up,'up',shape_up,hw_up,top_up,npl0,dgds0)
dn_ramp = ps.make_ramp(s_dn,'dn',shape_dn,hw_dn,top_dn,npl0,dgds0)
plasma  = ps.make_plasma(bulk,up_ramp,dn_ramp)

# define beam parameters
gbC    = 20000 # relativistic lorentz factor
eps    = 5e-6  # m-rad, normalized emittance
beta   = 0.10 # m
alpha  = 0.00
gamma  = (1.0+alpha**2)/beta # 1/m
dgb    = 0.01
dz     = 0
npart = 10
dist  = 'gauss'

# make beam
twiss0 = pb.make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz)
parts0 = pb.make_parts(twiss0,npart,dist)
ebeam  = pb.make_ebeam(0,twiss0,parts0)

# set beam waist position

# propagate beam backward from waist to start of plasma

# propagate beam through plasma

# analyze results

# plot results
