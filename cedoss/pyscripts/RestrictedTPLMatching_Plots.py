#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 09:42:42 2018

For restricted TPL matching, plots KL vs zm

@author: chris
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc

zm = np.linspace(-0.2, 0.2, 201)
bi = 0.02
bf = 0.06

kl = (2*bf*bi*zm + (bi+bf)*np.sqrt(bf*bi*(np.square(bi-bf)+np.square(zm)))) / (bf*bi*(np.square(bi+bf)+np.square(zm)))

#plt.plot(zm, kl)
#plt.grid(); plt.show()

gamma = Foc.gam_def
tpl_n = 0.5
tpl_l = Foc.Calc_Square_Lens(tpl_n*1e17, 1/kl*100, gamma)
    
plt.plot(zm*100, tpl_l, label=r'$\beta_i=$'+str(bi*100)+r'$\,\mathrm{[cm]}, \ \beta_f=$'+str(bf*100)+r'$\,\mathrm{[cm]}$')
plt.title("Required TPL thickness VS waist-waist separation")
plt.xlabel(r'$z_m \ \mathrm{[cm]}$')
plt.ylabel(r'$\mathrm{TPL \ Thickness \ [\mu m]}$')
plt.grid(); plt.legend(); plt.show()