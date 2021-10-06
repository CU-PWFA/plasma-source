#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 12:40:15 2021

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../")
from modules import ThreeDimensionAnalysis as ThrDim

np_0=np.array([2e16,1e16,2e16,2e16,2e16,2e16,3e17,3e17])#  ,2e16])
dndy=np.array([8e17,4e17,8e17,2e17,2.5e16,2.5e15,6e18,1.2e19])#  ,1])
rp_ave=np.array([80.83,103.7,80.83,80.83,80.83,80.83,37.4,37.4])
yoff=np.array([0.0409,0.0435,0.0403,0.00931,0.00124,7.98e-5,0.0117,0.0232])

rp_top=np.array([85.03,111.22,87.3,82.96,80.7,80.58,37.91,38.44])
rp_bot=np.array([76.51,96.92,78.35,80.74,80.48,80.55,36.91,36.47])
"""
facA = 0.3#0.33
#rptop_emp = np.sqrt((np_0+1/3*rp_ave*1e-4*dndy)/np_0)*rp_ave
rptop_emp = np.sqrt((np_0+facA*rp_ave*1e-4*dndy)/np_0)*rp_ave
rpbot_emp = np.sqrt((np_0-facA*rp_ave*1e-4*dndy)/np_0)*rp_ave
"""
#facB = 3.7#2.7             #2.7 1e16 and 2e16, 3.7 for 3e17
#kp = 1.883e-6 * np.sqrt(np_0) *1e-4
#yoff_emp = 1/2*(rptop_emp-rpbot_emp)/(facB*kp**-1)

lhalf = np.array([100,130,100,100,100,100,41,41])
yoff_emp = 1/2*(rptop_emp-rpbot_emp)/(lhalf)

#ns0 = np.array([7.07e16,7.08e16,3.78e16,6.96e16,7.10e16,2.00e18,8.50e17,7.08e16,7.08e16])
dndys = np.array([1.14e18,6.33e17,1.26e18,2.31e17,5.68e16,1.88e16,2.17e19,4.78e19])#  ,1])

#ns0t = np_0 * 3.5

#green set is 2e16
gset = np.where(np_0==2e16)[0]
yset = np.where(np_0==1e16)[0]
bset = np.where(np_0==3e17)[0]

"""
plt.title("R_top vs R_ave")
plt.xlabel(r'$R_{ave} \ (\mu m)$')
plt.ylabel(r'$R_{top} \ (\mu m)$')
plt.scatter(rp_ave[gset],rp_top[gset],c='g', label="2e16 cm-3")
plt.scatter(rp_ave[yset],rp_top[yset],c='y', label="1e16 cm-3")
plt.scatter(rp_ave[bset],rp_top[bset],c='b', label="3e17 cm-3")
plt.plot([30,110],[30,110],label="1:1", ls='--')
plt.grid(); plt.legend(); plt.show()

rt_fit = np.polyfit(rptop_emp, rp_top, 1)

plt.title("Emperical Guess for R_top")
plt.xlabel(r'$R_{top,emp} \ (\mu m)$')
plt.ylabel(r'$R_{top} \ (\mu m)$')
plt.scatter(rptop_emp[gset],rp_top[gset],c='g', label="2e16 cm-3")
plt.scatter(rptop_emp[yset],rp_top[yset],c='y', label="1e16 cm-3")
plt.scatter(rptop_emp[bset],rp_top[bset],c='b', label="3e17 cm-3")
plt.plot([30,110],[30,110],label="1:1", ls='--')
#plt.plot(rptop_emp, rptop_emp*rt_fit[0]+rt_fit[1], ls='--')
plt.grid(); plt.legend(); plt.show()
print("Fit slope = ",rt_fit[0])

rb_fit = np.polyfit(rpbot_emp, rp_bot, 1)

plt.title("Emperical Guess for R_bot")
plt.xlabel(r'$R_{bot,emp} \ (\mu m)$')
plt.ylabel(r'$R_{bot} \ (\mu m)$')
plt.scatter(rpbot_emp[gset],rp_bot[gset],c='g', label="2e16 cm-3")
plt.scatter(rpbot_emp[yset],rp_bot[yset],c='y', label="1e16 cm-3")
plt.scatter(rpbot_emp[bset],rp_bot[bset],c='b', label="3e17 cm-3")
plt.plot([30,110],[30,110],label="1:1", ls='--')
plt.grid(); plt.legend(); plt.show()
print("Fit slope = ",rb_fit[0])

plt.title("Error plot for Rtop and Rbot")
plt.xlabel("Simulation Case")
plt.ylabel("Error")
plt.plot(-(rp_top-rptop_emp)/rp_top,label='top')
plt.plot(-(rp_bot-rpbot_emp)/rp_bot,label='bot')
plt.grid(); plt.legend(); plt.show()



plt.title("Wave Offset Growth across Simulations")
plt.xlabel(r'$\partial n / \partial y \ (cm^{-4})$')
plt.ylabel(r'$y_{off} \, / \, \mu m$')
plt.scatter(dndy[gset],yoff[gset],c='g', label="2e16 cm-3")
plt.scatter(dndy[yset],yoff[yset],c='y', label="1e16 cm-3")
plt.scatter(dndy[bset],yoff[bset],c='b', label="3e17 cm-3")
#plt.plot([30,110],[30,110],label="1:1", ls='--')
plt.grid(); plt.legend(); plt.show()

yoff_fit = np.polyfit(yoff_emp, yoff, 1)
plt.title("Emperical Guess for y_off")
plt.xlabel(r'$y_{off,emp} \, / \, \mu m$')
plt.ylabel(r'$y_{off} \, / \, \mu m$')
plt.scatter(yoff_emp[gset],yoff[gset],c='g', label="2e16 cm-3")
plt.scatter(yoff_emp[yset],yoff[yset],c='y', label="1e16 cm-3")
plt.scatter(yoff_emp[bset],yoff[bset],c='b', label="3e17 cm-3")
plt.plot([-0.005,0.05],[-0.005,0.05],label="1:1", ls='--')
plt.grid(); plt.legend(); plt.show()
print("Fit slope = ",yoff_fit[0])

nz = np.where(yoff!=0)[0]
zs = np.where(yoff==0)[0]
plt.title("Error plot for nonzero yOff")
plt.xlabel("Simulation Case")
plt.ylabel("Error")
plt.plot(-(yoff[nz]-yoff_emp[nz])/yoff[nz])
plt.grid(); plt.show()
"""


x = np.logspace(np.log10(min(dndy)),np.log10(max(dndy)),50)
power = 3
dndys_fit = np.polyfit(dndy, dndys, power)
"""
dndys_fit[0]=-3.0789e-38
dndys_fit[1]=6.1568e-19
dndys_fit[2]=1.0266
dndys_fit[3]=2.5089e16
"""
y = np.zeros(len(x))
print("Fit with a power of",power,":")
print("(dn/dy)_sheath = ")
for i in range(power+1):
    y = y+ np.power(x,power-i)*dndys_fit[i]
    print("    (dn/dy)^",power-i, " x " , dndys_fit[i]," +")

plt.title("Sheath Density Gradient")
plt.xlabel(r'$\partial n / \partial y \ (cm^{-4})$')
plt.ylabel(r'$\partial n / \partial y_s \ (cm^{-4})$')
plt.loglog(x,y,ls='--',label="Third Order Fit")
plt.scatter(dndy[gset],dndys[gset],c='g', label="2e16 cm-3")
plt.scatter(dndy[yset],dndys[yset],c='y', label="1e16 cm-3")
plt.scatter(dndy[bset],dndys[bset],c='b', label="3e17 cm-3")
#plt.plot([30,110],[30,110],label="1:1", ls='--')
plt.grid(); plt.legend(); plt.show()

"""#Don't care, doesn't show up, has longitudinal variation, idk
ns0_fit = np.polyfit(np_0, ns0, 1)
plt.title("Sheath Central Density")
plt.xlabel(r'$ n_p \times 3.5 \ (cm^{-3})$')
plt.ylabel(r'$ n_s \ (cm^{-3})$')
plt.loglog((min(ns0),max(ns0)),[min(ns0),max(ns0)],label="1:1",ls='--')
plt.scatter(ns0t[gset],ns0[gset],c='g', label="2e16 cm-3")
plt.scatter(ns0t[yset],ns0[yset],c='y', label="1e16 cm-3")
plt.scatter(ns0t[bset],ns0[bset],c='b', label="3e17 cm-3")
#plt.plot([30,110],[30,110],label="1:1", ls='--')
plt.grid(); plt.legend(); plt.show()
"""
"""
yoff_fit = np.polyfit(yoff_emp, yoff, 1)
plt.title("Emperical Guess for y_off")
plt.xlabel(r'$y_{off,emp} \, / \, \mu m$')
plt.ylabel(r'$y_{off} \, / \, \mu m$')
plt.scatter(yoff_emp[gset],yoff[gset],c='g', label="2e16 cm-3")
plt.scatter(yoff_emp[yset],yoff[yset],c='y', label="1e16 cm-3")
plt.scatter(yoff_emp[bset],yoff[bset],c='b', label="3e17 cm-3")
plt.plot([-0.005,0.05],[-0.005,0.05],label="1:1", ls='--')
plt.grid(); plt.legend(); plt.show()
print("Fit slope = ",yoff_fit[0])

nz = np.where(yoff!=0)[0]
zs = np.where(yoff==0)[0]
plt.title("Error plot for nonzero yOff")
plt.xlabel("Simulation Case")
plt.ylabel("Error")
plt.plot(-(yoff[nz]-yoff_emp[nz])/yoff[nz])
plt.grid(); plt.show()
"""