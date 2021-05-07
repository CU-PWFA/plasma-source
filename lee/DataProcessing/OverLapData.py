#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 18:31:51 2020

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt
#%%
xsolved_14= (np.load('Data14_R.npy')[1]-np.load('Data14_L.npy')[1])/(np.load('Data14_L.npy')[0]-np.load('Data14_R.npy')[0])
ysolved_14= np.load('Data14_L.npy')[0]*xsolved_14+np.load('Data14_L.npy')[1]

xsolved_20= (np.load('Data20_R.npy')[1]-np.load('Data20_L.npy')[1])/(np.load('Data20_L.npy')[0]-np.load('Data20_R.npy')[0])
ysolved_20= np.load('Data20_L.npy')[0]*xsolved_20+np.load('Data20_L.npy')[1]

xsolved_24= (np.load('Data24_R.npy')[1]-np.load('Data24_L.npy')[1])/(np.load('Data24_L.npy')[0]-np.load('Data24_R.npy')[0])
ysolved_24= np.load('Data24_L.npy')[0]*xsolved_24+np.load('Data24_L.npy')[1]

xsolved_26= (np.load('Data26_R.npy')[1]-np.load('Data26_L.npy')[1])/(np.load('Data26_L.npy')[0]-np.load('Data26_R.npy')[0])
ysolved_26= np.load('Data26_L.npy')[0]*xsolved_26+np.load('Data26_L.npy')[1]

xsolved_28= (np.load('Data28_R.npy')[1]-np.load('Data28_L.npy')[1])/(np.load('Data28_L.npy')[0]-np.load('Data28_R.npy')[0])
ysolved_28= np.load('Data28_L.npy')[0]*xsolved_28+np.load('Data28_L.npy')[1]

xsolved_32= (np.load('Data32_R.npy')[1]-np.load('Data32_L.npy')[1])/(np.load('Data32_L.npy')[0]-np.load('Data32_R.npy')[0])
ysolved_32= np.load('Data32_L.npy')[0]*xsolved_32+np.load('Data32_L.npy')[1]

xsolved_38= (np.load('Data38_R.npy')[1]-np.load('Data38_L.npy')[1])/(np.load('Data38_L.npy')[0]-np.load('Data38_R.npy')[0])
ysolved_38= np.load('Data38_L.npy')[0]*xsolved_38+np.load('Data38_L.npy')[1]

xsolved_42= (np.load('Data42_R.npy')[1]-np.load('Data42_L.npy')[1])/(np.load('Data42_L.npy')[0]-np.load('Data42_R.npy')[0])
ysolved_42= np.load('Data42_L.npy')[0]*xsolved_42+np.load('Data42_L.npy')[1]

xsolved_48= (np.load('Data48_R.npy')[1]-np.load('Data48_L.npy')[1])/(np.load('Data48_L.npy')[0]-np.load('Data48_R.npy')[0])
ysolved_48= np.load('Data48_L.npy')[0]*xsolved_48+np.load('Data48_L.npy')[1]
#%%
x= np.linspace(-1000, 2000, 2000)
Data14_L= np.load('Data14_L.npy')[0]*(x+xsolved_14)+np.load('Data14_L.npy')[1]-ysolved_14
Data14_R= np.load('Data14_R.npy')[0]*(x+xsolved_14)+np.load('Data14_R.npy')[1]-ysolved_14

Data20_L= np.load('Data20_L.npy')[0]*(x+xsolved_20)+np.load('Data20_L.npy')[1]-ysolved_20
Data20_R= np.load('Data20_R.npy')[0]*(x+xsolved_20)+np.load('Data20_R.npy')[1]-ysolved_20

Data24_L= np.load('Data24_L.npy')[0]*(x+xsolved_24)+np.load('Data24_L.npy')[1]-ysolved_24
Data24_R= np.load('Data24_R.npy')[0]*(x+xsolved_24)+np.load('Data24_R.npy')[1]-ysolved_24

Data26_L= np.load('Data26_L.npy')[0]*(x+xsolved_26)+np.load('Data26_L.npy')[1]-ysolved_26
Data26_R= np.load('Data26_R.npy')[0]*(x+xsolved_26)+np.load('Data26_R.npy')[1]-ysolved_26

Data28_L= np.load('Data28_L.npy')[0]*(x+xsolved_28)+np.load('Data28_L.npy')[1]-ysolved_28
Data28_R= np.load('Data28_R.npy')[0]*(x+xsolved_28)+np.load('Data28_R.npy')[1]-ysolved_28

Data32_L= np.load('Data32_L.npy')[0]*(x+xsolved_32)+np.load('Data32_L.npy')[1]-ysolved_32
Data32_R= np.load('Data32_R.npy')[0]*(x+xsolved_32)+np.load('Data32_R.npy')[1]-ysolved_32

Data38_L= np.load('Data38_L.npy')[0]*(x+xsolved_38)+np.load('Data38_L.npy')[1]-ysolved_38
Data38_R= np.load('Data38_R.npy')[0]*(x+xsolved_38)+np.load('Data38_R.npy')[1]-ysolved_38

Data42_L= np.load('Data42_L.npy')[0]*(x+xsolved_42)+np.load('Data42_L.npy')[1]-ysolved_42
Data42_R= np.load('Data42_R.npy')[0]*(x+xsolved_42)+np.load('Data42_R.npy')[1]-ysolved_42

Data48_L= np.load('Data48_L.npy')[0]*(x+xsolved_48)+np.load('Data48_L.npy')[1]-ysolved_48
Data48_R= np.load('Data48_R.npy')[0]*(x+xsolved_48)+np.load('Data48_R.npy')[1]-ysolved_48
#%%
theta_14= np.arctan(abs(Data14_R[0]/abs(x)[0]))
Dy_14L= x*np.sin(theta_14)+(np.load('Data14_L.npy')[0]*x)*np.cos(theta_14)
Dy_14R= x*np.sin(theta_14)+(np.load('Data14_R.npy')[0]*x)*np.cos(theta_14)

theta_20= np.arctan(abs(Data20_R[0]/abs(x)[0]))
Dy_20L= x*np.sin(theta_20)+(np.load('Data20_L.npy')[0]*x)*np.cos(theta_20)
Dy_20R= x*np.sin(theta_20)+(np.load('Data20_R.npy')[0]*x)*np.cos(theta_20)

theta_24= np.arctan(abs(Data24_R[0]/abs(x)[0]))
Dy_24L= x*np.sin(theta_24)+(np.load('Data24_L.npy')[0]*x)*np.cos(theta_24)
Dy_24R= x*np.sin(theta_24)+(np.load('Data24_R.npy')[0]*x)*np.cos(theta_24)

theta_26= np.arctan(abs(Data26_R[0]/abs(x)[0]))
Dy_26L= x*np.sin(theta_26)+(np.load('Data26_L.npy')[0]*x)*np.cos(theta_26)
Dy_26R= x*np.sin(theta_26)+(np.load('Data26_R.npy')[0]*x)*np.cos(theta_26)

theta_28= np.arctan(abs(Data28_R[0]/abs(x)[0]))
Dy_28L= x*np.sin(theta_28)+(np.load('Data28_L.npy')[0]*x)*np.cos(theta_28)
Dy_28R= x*np.sin(theta_28)+(np.load('Data28_R.npy')[0]*x)*np.cos(theta_28)

theta_32= np.arctan(abs(Data32_R[0]/abs(x)[0]))
Dy_32L= x*np.sin(theta_32)+(np.load('Data32_L.npy')[0]*x)*np.cos(theta_32)
Dy_32R= x*np.sin(theta_32)+(np.load('Data32_R.npy')[0]*x)*np.cos(theta_32)

theta_38= np.arctan(abs(Data38_R[0]/abs(x)[0]))
Dy_38L= x*np.sin(theta_38)+(np.load('Data38_L.npy')[0]*x)*np.cos(theta_38)
Dy_38R= x*np.sin(theta_38)+(np.load('Data38_R.npy')[0]*x)*np.cos(theta_38)

theta_42= np.arctan(abs(Data42_R[0]/abs(x)[0]))
Dy_42L= x*np.sin(theta_42)+(np.load('Data42_L.npy')[0]*x)*np.cos(theta_42)
Dy_42R= x*np.sin(theta_42)+(np.load('Data42_R.npy')[0]*x)*np.cos(theta_42)

theta_48= np.arctan(abs(Data48_R[0]/abs(x)[0]))
Dy_48L= x*np.sin(theta_48)+(np.load('Data48_L.npy')[0]*x)*np.cos(theta_48)
Dy_48R= x*np.sin(theta_48)+(np.load('Data48_R.npy')[0]*x)*np.cos(theta_48)

#%%
tan_14= abs(Dy_14L[0])/abs(x[0])
tan_20= abs(Dy_20L[0])/abs(x[0])
tan_24= abs(Dy_24L[0])/abs(x[0])
tan_26= abs(Dy_26L[0])/abs(x[0])
tan_28= abs(Dy_28L[0])/abs(x[0])
tan_32= abs(Dy_32L[0])/abs(x[0])
tan_38= abs(Dy_38L[0])/abs(x[0])
tan_42= abs(Dy_42L[0])/abs(x[0])
tan_48= abs(Dy_48L[0])/abs(x[0])
#%%
plt.plot(x, Dy_14R, 'k')
plt.plot(x, Dy_14L, 'b', label='6cm')
plt.plot(x, Dy_20L, 'g', label= '12cm')
plt.plot(x, Dy_24L, 'silver', label= '16m')
plt.plot(x, Dy_26L, 'r', label='31.8cm' )
plt.plot(x, Dy_28L, 'limegreen', label= '34.8cm')
plt.plot(x, Dy_32L, 'c', label='40.8cm')
plt.plot(x, Dy_38L, 'm', label='49.8cm')
plt.plot(x, Dy_42L, 'y', label='71.1cm')
#plt.plot(x, Dy_48L, 'darkorange', label='86cm')
plt.xlabel('x (pixel)')
plt.ylabel('y (pixel)')
#plt.axis('scaled')
plt.legend()
#%%
plt.plot(6, tan_14, '.')
plt.plot(12, tan_20, '.')
plt.plot(16, tan_24, '.')
plt.plot(31.8, tan_26, '.')
plt.plot(34.8, tan_28, '.')
plt.plot(40.8, tan_32, '.')
plt.plot(49.8, tan_38, '.')
plt.plot(71.1, tan_42, '.')
#plt.plot(86, tan_48, '.')

#plt.xlabel('x (pixel)')
#plt.ylabel('y (pixel)')
#plt.axis('scaled')
#plt.legend()
#%%
plt.plot(x, Data28_L)
plt.plot(x, Data28_R)
plt.plot(x, Data24_L)
plt.plot(x, Data24_R)
