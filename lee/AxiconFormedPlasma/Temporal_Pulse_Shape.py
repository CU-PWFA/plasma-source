# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 13:22:51 2020

@author: Robert
"""
#%%
ind = 99
Nt = pulse.Nt
T = pulse.T
Nx = pulse.Nx
X = pulse.X
e = np.zeros((Nt, Nx), dtype='complex128')
e[:, :] = pulse.load_field(ind)[0]
I = ionization.intensity_from_field(e)
I_max = np.amax(I)
dx = X/(Nx-1)

ext = [-T/2, T/2, -X/2, X/2-dx]
plt.figure(figsize=(4, 2), dpi=300)
plt.imshow(np.fliplr(np.flipud(np.transpose(I))), aspect='auto', extent=ext, cmap='viridis', 
           interpolation='Spline16')
cb = plt.colorbar(format="%0.2f")
cb.set_label(r'Laser Intensity ($10^{14} W/cm^2$)')
plt.xlabel('t (fs)')
plt.ylabel(r'x ($\mathrm{\mu m}$)')
plt.ylim([-500, 500])
plt.show()