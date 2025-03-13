#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 11:21:40 2023

Stating some of the beam paramters from the Aug E308 2022 Shift

@author: chris
"""

import numpy as np

#Sigmas from before and after, averaged
sigx = 38.55
sigx_err = 1.69
sigy = 32.50
sigy_err = 2.25

#Sig_z's from TCAV.  Before and After, averaged
sigz = 27.86
sigz_err = 5.12

#Energy from multi-wire scan (?) Before and after, average
energy = 9.987
energy_approx = 10

#Emittances from multi-wire scans.  Only before.  The after is a lil late?
emitxn = 23.63
emitxn_err = 18.23
emityn = 10.46
emityn_err = 2.86

#Charge, in nC
charge = 1.5

#Designed betafunctions at focus
betax_des = 0.06
betay_des = 0.06

#Calculating the peak density of the beam, in cm^-3
e=1.6022e-19
rho_peak = 1/e*charge*1e-9/(np.power(2*np.pi,3/2)*sigx*sigy*sigz*np.power(1e-6,3))/np.power(100,3)
print("Peak Density: ",rho_peak)

#error?
factors = 1/e*charge*1e-9/(np.power(2*np.pi,3/2)*np.power(1e-6,3))/np.power(100,3)
error = np.sqrt(factors**2*(sigx_err**2/sigx**4/sigy**2/sigz**2
                            +sigy_err**2/sigy**4/sigx**2/sigz**2
                            +sigz_err**2/sigz**4/sigx**2/sigy**2))
print("Error:",error);print()

gam = 19569.5
betax = (sigx*1e-6)**2*gam/(emitxn*1e-6)
betay = (sigy*1e-6)**2*gam/(emityn*1e-6)

betax_err = np.sqrt((2*sigx*1e-6*gam/(emitxn*1e-6)*sigx_err*1e-6)**2
                    +((sigx*1e-6)**2*gam/(emitxn*1e-6)**2*emitxn_err*1e-6)**2)
betay_err = np.sqrt((2*sigy*1e-6*gam/(emityn*1e-6)*sigy_err*1e-6)**2
                    +((sigy*1e-6)**2*gam/(emityn*1e-6)**2*emityn_err*1e-6)**2)
print("Calculated Betax:",betax)
print(" Error: ",betax_err)
print("Calculated Betay:",betay)
print(" Error: ",betay_err)



