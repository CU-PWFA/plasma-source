#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:13:08 2020

Looking at some phase space plots of the "slingshotted" plasma electrons in the
last dump of plasma lens PIC simulations

@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np
import timeit
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "../../python")
from vsim import load

path = '/home/chris/Desktop/BeamProp/vorpaltest'
runpath = '/home/chris/Desktop/200e16/'

filename = runpath + 'PTPLDoubleTanh_electrons_10.h5'
rhofile = runpath + 'PTPLDoubleTanh_rhoPlasma_10.h5'

rho = np.sum(np.array(load.get_field_data(rhofile,'rhoPlasma')))
dx = 5.629629629629629e-07; dy = 5.639097744360901e-07; dz = 5.639097744360901e-07;
total_charge = rho * dx * dy * dz
print("Charge = ",total_charge*1e9,"nC")

debug = 0

dump = 20
cores = 4
threshold = 0.000001

beam = PProp.VorpalBeam(path, filename, threshold, debug=debug)
print("N: ",beam.N)

mingm = 1.2
maxrad = 1000 *1e-3 #mrad
minz = -1000

ptcls = beam.load_ptcls(0)[0]
weights = ptcls[:,6]
presum = np.sum(weights)

sort = np.argsort(weights)
ptcls = ptcls[sort]

gmcutoff = np.array(np.where(ptcls[:,5] >= mingm)[0])
ptcls = ptcls[gmcutoff]

radcutoffx = np.array(np.where(np.abs(ptcls[:,3]) <= maxrad)[0])
ptcls = ptcls[radcutoffx]
radcutoffy = np.array(np.where(np.abs(ptcls[:,1]) <= maxrad)[0])
ptcls = ptcls[radcutoffy]

minzcuttof = np.array(np.where(ptcls[:,4]*1e6 >= minz)[0])
ptcls = ptcls[minzcuttof]

ptcls_x = ptcls[:,2]*1e6
ptcls_xp= ptcls[:,3]*1e3
ptcls_y = ptcls[:,0]*1e6
ptcls_yp= ptcls[:,1]*1e3
ptcls_z = ptcls[:,4]*1e6
ptcls_gm= ptcls[:,5]
weights = ptcls[:,6]

print("Average gamma:" ,np.average(ptcls_gm,weights=weights))
print("Max gamma:" ,np.max(ptcls_gm))

postsum = np.sum(weights)
postcharge = total_charge*postsum/presum
print("Plotted Charge = ",postcharge*1e9,"nC")

numbins = 101   #bins in 1d hist
binno = 120      #bins in 2d hist
xrange = 10#8
yrange = 80.0

#plt.hist(ptcls_gm,weights=weights,bins=numbins,log=True)
#plt.show()

plt.title("z-"+r'$\gamma$'+" plane")
plt.hist2d(ptcls_z, ptcls_gm, weights=weights, bins=(binno,binno), cmap=plt.cm.jet)#, range = [[-xrange,xrange], [-yrange,yrange]])
plt.xlabel("z"+r'$\mathrm{\ [\mu m]}$')
plt.ylabel(r'$\gamma$')
plt.show()

#sys.exit()

plt.title("z-x plane")
plt.hist2d(ptcls_z, ptcls_x, weights=weights, bins=(binno,binno), cmap=plt.cm.jet)#, range = [[-xrange,xrange], [-yrange,yrange]])
plt.xlabel("z"+r'$\mathrm{\ [\mu m]}$')
plt.ylabel("x"+r'$\mathrm{\ [\mu m]}$')
plt.show()

plt.title("z-y plane")
plt.hist2d(ptcls_z, ptcls_y, weights=weights, bins=(binno,binno), cmap=plt.cm.jet)#, range = [[-xrange,xrange], [-yrange,yrange]])
plt.xlabel("z"+r'$\mathrm{\ [\mu m]}$')
plt.ylabel("y"+r'$\mathrm{\ [\mu m]}$')
plt.show()

plt.title("x-y plane")
plt.hist2d(ptcls_x, ptcls_y, weights=weights, bins=(binno,binno), cmap=plt.cm.jet)#, range = [[-xrange,xrange], [-yrange,yrange]])
plt.xlabel("x"+r'$\mathrm{\ [\mu m]}$')
plt.ylabel("y"+r'$\mathrm{\ [\mu m]}$')
plt.show()

plt.title("x'-y' plane")
plt.hist2d(ptcls_xp, ptcls_yp, weights=weights, bins=(binno,binno), cmap=plt.cm.jet)#, range = [[-xrange,xrange], [-yrange,yrange]])
plt.xlabel("x'"+r'$\mathrm{\ [mm-mrad]}$')
plt.ylabel("y'"+r'$\mathrm{\ [mm-mrad]}$')
plt.show()

plt.title("Quick Template")
plt.hist2d(ptcls_gm, ptcls_y, weights=weights, bins=(binno,binno), cmap=plt.cm.jet)#, range = [[-xrange,xrange], [-yrange,yrange]])
plt.show()

#beam.plot_phase_hist_at(0, fitted=True)
#print(beam.get_emittance_n(0)[0],beam.get_emittance_n(0)[1])
#print("Full emittance: ",np.average(beam.get_emittance_n(0)))