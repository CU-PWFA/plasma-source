#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 19:30:26 2021

@author: valentinalee
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 17:27:03 2020

@author: valentinalee
"""

#%%'

import sys
import os
sys.path.insert(0, "python")
import pyximport; pyximport.install()
import numpy as np
from beam.beams import laserbeam
from beam.elements import plasma
from beam import interactions
from ionization import ionization
from lens import profile
from lens import phaselens
from lens import ray
import matplotlib.pyplot as plt
from unwrap import unwrap
from propagation import laser
from scipy import interpolate;
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
'''
If laser plasma interaction fail, check the grid density. If the spacing k*theta is larger than the grid spacig, it will fail.
Try with a smaller crossing angle to confirm the code works. 
'''

os.chdir("../")
#%%
c= 3e8
epsilon= 8.85e-12

laserPara= np.load('laser_para_88.npy')
FFphase=np.load('FFphase_88.npy')
PhaseGrid=np.load('grid_88.npy')
f= interp2d(PhaseGrid*1e6, PhaseGrid*1e6, FFphase, kind='cubic')

#%%
PlasmaWidthS=45
density= 1

basePath = 'lee/'
path = basePath + 'LaserRefraction/'

power= 5
# Build the plasma
plasmaParams = {
    'Nx' : 2**14, #2**16, #2**17,
    'Ny' : 2**8, #2**8,
    'Nz' : 2**10, #2**11,
    'X' : 15000,#0.060e6,
    'Y' : 3200,#6000, <- it works when 6000
    'Z' : 0.14e6, 
    'n0' : density, #Nominal gas number density in 10^17 cm^-3.
    'atom' : ionization.He,
    'path' : path,
    'name' : 'HePlasma',
    'load' : False,
    'cyl' : True
}

w = PlasmaWidthS
z0 = -0.08 #The distance at which the uniform fully ionized plasma starts.
zf = plasmaParams['Z']
dz = 0.3e6#0.3e6 #The length of the fully ionized plasma.
sigmaIn = 0.05e6#0.15e6
sigmaOut = 0.05e6#0.15e6
z, nez = profile.plasma_gaussian_ramps(z0, dz, sigmaIn, sigmaOut, plasmaParams['Nz'], zf)
He = plasma.Plasma(plasmaParams)
r = np.sqrt(He.x[:, None]**2 + He.y[None, :]**2)

#%
os.chdir("lee/LaserRefraction/elements/element_HePlasma")
for z in range(plasmaParams['Nz']):
    print(z)
    ne = plasmaParams['n0'] * nez[z] * np.exp(-(r**2/(2*w**2))**power)
    ne= ne[:, int(plasmaParams['Ny']/2)]
    np.save(str(plasmaParams['name']) + '_plasmaDensity_' + str(z) + '.npy', ne)

os.chdir("../../../../")
    
del nez
del ne
#%%
He.plot_long_density_center()
#%%
beamParams = {
    'Nx' : plasmaParams['Nx'],
    'Ny' : plasmaParams['Ny'],
    'X' : plasmaParams['X'],
    'Y' : plasmaParams['Y'],
    'lam' : 0.8, #wavelength
    'path' : path,
    'name' : 'ProbePulse',
    'load' : False,
    'threads' : 8,
    'cyl' : False,
    'E0' : laserPara[0], #5.473494, # Gives peak input intensity of 3.98e10, that is, E=0.2mJ, beam dia= 8mm, pulse width= 100ps
    'waist' : laserPara[1]*1e6,
    'z0' : 0.0,
    'order' : laserPara[2],
    'theta' : 4, #3#The angle of propagation in degrees, positive is angled upward.
    'dx' : -4.88e3
}

beam = laserbeam.GeneralSuperGaussianLaser(beamParams)
#%
#xgrid= np.linspace(-beamParams['X']/2, beamParams['X']/2, beamParams['Nx'], dtype='double')
xgrid= np.linspace(-beamParams['X']/2-beamParams['dx'], beamParams['X']/2-beamParams['dx'], beamParams['Nx'], dtype='double')
ygrid= np.linspace(-beamParams['Y']/2, beamParams['Y']/2, beamParams['Ny'], dtype='double')

phase2= f(xgrid, ygrid)

#%
beam.set_field(beam.e * np.exp(-1j*np.transpose(phase2)))
#%
#interactions.beam_phase(beam, phase2)
#%
I= abs(beam.e)**2*c*epsilon/2

#%
#xgridP= np.linspace(-beamParams['X']/2, beamParams['X']/2, beamParams['Nx'], dtype='double')
#%
#propE= laser.fourier_prop2(beam.e, xgridP, ygrid, 0.8e6, 0.8)
#%
#plt.figure(15)
#plt.axis('scaled')
#plt.pcolormesh(xgridP, ygrid, np.transpose(abs(beam.e)**2))
#plt.pcolormesh(xgridP, ygrid, np.transpose(abs(propE[0, :, :])**2))
#plt.pcolormesh(xgridP, ygrid, (phase2))
#%
interactions.beam_plasma(beam, He)
#%%
plt.figure(2)   
beam.plot_intensity_at(1023)
#plt.savefig('Int')
#plt.close()

#%%           
#Fix propagation angle
NewBeam = beam.e#*np.exp(-1j*beam.k*np.radians(beam.theta)*beam.x[:, None])
#NewBeam = beam.e*np.exp(-1j*beam.k*np.radians(beam.theta)*beam.x[:, None])
x_final= plasmaParams['Z']*np.tan(np.radians(beamParams['theta']))+beamParams['dx']
xpixel= beamParams['X']/beamParams['Nx']
ypixel= beamParams['Y']/beamParams['Ny']

nth_pixel= round((x_final/xpixel)+beamParams['Nx']/2)

indexSY=0
indexFY=beamParams['Ny']
indexSX=int(nth_pixel-abs(2*(round(beamParams['waist']/xpixel))))
indexFX=int(nth_pixel+abs(2*(round(beamParams['waist']/xpixel))))
E_afterPlasma= np.swapaxes(np.asarray(NewBeam[indexSX:indexFX, indexSY:indexFY]), 0, 1)

#%%
print('PW'+str(PlasmaWidthS)+'_GP'+str(power)+'_n'+str(density))


np.save('PW'+str(PlasmaWidthS)+'_GP'+str(power)+'_n'+str(int(density*10)), E_afterPlasma)
            
