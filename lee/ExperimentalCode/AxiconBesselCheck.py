#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 11:44:25 2021

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import scipy
#%%
Axicon= np.array(Image.open('19423598_2109240001_0004.tiff'))[325:575, 1780:1980]

#%%
plt.pcolormesh(Axicon)

#%%
x= np.linspace(0, 200, 1000)
bessel= scipy.special.jv(0, (x-100)*2)
plt.plot(x, (bessel)**2)

#%%
#plt.plot(Axicon[124, :])
plt.plot(Axicon[:, 97])
#%%                
xvar= np.linspace(0, 200, 200)
yvar= Axicon[124, :]/np.amax(Axicon[124, :])
fitfn= lambda p, xvar: scipy.special.jv(0, (xvar+p[0])*p[1])**2

errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar
p0= [-100, 0.2]
p1, success= scipy.optimize.leastsq(errfunc, p0[:], args=(xvar, yvar), epsfcn= 1e-50)
print(p1)
Newx= np.linspace(0, 200, 1000)
#%%
plt.plot(xvar, yvar, '.')
plt.plot(Newx, fitfn(p1, Newx))

#%%                
xvar= np.linspace(0, 250, 250)
yvar= Axicon[:, 97]/np.amax(Axicon[:, 97])
fitfn= lambda p, xvar: scipy.special.jv(0, (xvar+p[0])*p[1])**2

errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar
p0= [-125, 0.2]
p1, success= scipy.optimize.leastsq(errfunc, p0[:], args=(xvar, yvar), epsfcn= 1e-50)
print(p1)
Newx= np.linspace(0, 200, 1000)
#%%
plt.plot(xvar, yvar, '.')
plt.plot(Newx, fitfn(p1, Newx))
