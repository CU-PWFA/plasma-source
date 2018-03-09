#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 19:09:17 2018

@author: keenan
"""

import KH_PlasmaProp as PProp
from collections import defaultdict
from matplotlib import pyplot as plt
import numpy as np
import beam_ana as ba
import particle_beam_propagation as pbp
import particle_beam as pb
def PlotExit(exit_beam, plasma):
    nstep      = len(exit_beam)
    s_exit     = np.zeros(nstep)
    exit_sig   = np.zeros(nstep)
    exit_x_eps = np.zeros(nstep)
    exit_beta  = np.zeros(nstep)
    frac       = 1.00
    for i in range(0,nstep):
        exit_rms     = ba.calc_ebeam_rms(exit_beam, i, frac)
        s_exit[i]    = exit_beam[i]["s"]
        exit_beta[i] = exit_beam[i]["beta"]
        exit_x_eps[i] = exit_rms["x_eps"]
        exit_gbC     = exit_beam[i]["gbC"]
        exit_sig[i]  = np.sqrt(exit_x_eps[i]*exit_beta[i]/exit_gbC)*1e6
    figA, (ax1, ax3) = plt.subplots(2, sharex=True, sharey=False)
    # Plot beta for virtual and real beam
    ax1.plot(s_exit, exit_sig, color = 'b', linestyle = '-.')
    ax1.set_xlim([0.5,3.0])
    ax1.set_ylim([0,50])
    ax1.set_ylabel(r'$\sigma$ [$\mu$m]',color='b')
    ax1.tick_params('y',colors='b')   
   
params = PProp.ReturnDefaultParams()


params['npart'] = 100 #Want 10,000 to 100,000 for a final figure

twiss   = PProp.CallMakeTwiss(params)
parts   = PProp.CallMakeParts(twiss, params)
ebeam0  = PProp.CallMakeBeam(twiss, parts, params)
initialbeam = ebeam0.copy()
ebeam0  = PProp.PropagateBackwards(ebeam0, params)
plasma0 = PProp.MakeBulkPlasma(params)
vbeam0  = ebeam0.copy()
ebeam0  = PProp.PropagatePlasma(ebeam0, plasma0)
vbeam0  = PProp.PropagateVirtual(vbeam0, plasma0)
#Get virtual waist at exit
exitbeam           = defaultdict(dict)
exitbeam[0]        = ebeam0[len(ebeam0)-1]
exitbeam[0]['dgb'] = 0; exitbeam[0]['dz'] = 0
exit_waist            = 3.5 - (-exitbeam[0]['alpha'] / exitbeam[0]['gamma'])
ds                 = ebeam0[1]['s'] - ebeam0[0]['s'];
s                  = [3.5, 0]
pbp.prop_ebeam_drift(exitbeam,s)    
twiss = pb.get_twiss(exitbeam, len(exitbeam)-1);
parts = pb.get_parts(exitbeam,len(exitbeam)-1)
exitbeam = pb.make_ebeam(0,twiss[len(exitbeam)-1],parts[len(exitbeam)-1])
#plasma_s = np.arange(s_waist, 3.5, ds, dtype = 'float')
#exit_plasma = defaultdict(dict)
exitbeam    = PProp.PropagateVirtual(exitbeam, plasma0)
PProp.PlotPropagation(ebeam0, vbeam0, exitbeam, plasma0, log = True, waist = False)
#PlotExit(exitbeam, plasma0)