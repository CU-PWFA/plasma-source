#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 13:42:13 2021

Using the August 2021 runs

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../")
from modules import ThreeDimensionAnalysis as ThrDim
import numpy.polynomial.polynomial as poly

"""
#This set includes the control cases
np_0=np.array([2e16,2e16,2e16,2e16,2e16,1e17,1e17,1e17,1e16,1e16])
dndy=np.array([0,8e17,2e17,2.5e16,2.5e15,0,4e18,2e18,0,4e17])
rp_ave=np.array([76.15,76.15,76.15,76.15,76.15,48.67,48.67,48.67,91.86,91.86])
yoff=np.array([0,0.0346,0.00866,0.00107,0.000129,0,0.0261,0.0121,0,0.0403])

rp_top=np.array([76.15,80.41,77.09,76.24,76.171,48.67,50.56,49.59,91.86,98.04])
rp_bot=np.array([76.15,72.58,75.20,76.07,76.137,48.67,47.03,47.84,91.86,86.97])
lhalf = np.array([114,114,114,114,114,61.5,61.5,61.5,131.7,131.7])
dndys = np.array([0,8.64e17,1.47e17,6.01e16,3.96e16,0,1.21e19,6.54e18,0,7.96e17])
ering = np.array([0,8.833e8,3.0025e8,1.2321e8,1.0337e8,0,1.8627e9,1.1623e9,0,5.8867e8])
"""

#This set does NOT include the control cases
np_0=np.array([2e16,2e16,2e16,2e16,1e17,1e17,1e16])
dndy=np.array([8e17,2e17,2.5e16,2.5e15,4e18,2e18,4e17])
rp_ave=np.array([76.15,76.15,76.15,76.15,48.67,48.67,91.86])
yoff=np.array([0.0346,0.00866,0.00107,0.000129,0.0261,0.0121,0.0403])
garr = np.array([0.305,0.076,0.010,0.001,0.065,0.032,0.367])

rp_top=np.array([80.41,77.09,76.24,76.171,50.56,49.59,98.04])
rp_bot=np.array([72.58,75.20,76.07,76.137,47.03,47.84,86.97])
lhalf = np.array([114,114,114,114,61.5,61.5,131.7])
#FROM BEFORE FIX dndys = np.array([8.64e17,1.47e17,6.01e16,3.96e16,1.21e19,6.54e18,7.96e17])
#FROM BEFORE FIX dndys_std = np.array([1.90e17,1.34e17,5.43e16,2.65e16,4.93e18,2.81e18,1.88e17])
dndys = np.array([1.76e18,4.33e17,6.05e16,6.07e15,1.12e19,4.96e18,7.96e17])
dndys_std = np.array([1.18e17,4.44e16,4.02e16,4.65e15,3.03e18,3.15e18,1.88e17])
#ering = np.array([8.833e8,3.0025e8,1.2321e8,1.0337e8,1.8627e9,1.1623e9,5.8867e8])
#ering = np.array([7.808e8,1.867e8,2.048e7,9.196e5,1.608e9,9.568e8,5.429e8])

ering_norm = np.array([1.33e17,3.52e16,4.18e15,4.44e14,7.12e17,3.69e17,5.98e16])
ering_norm_std = np.array([7.86e15,3.91e15,1.64e15,4.28e14,5.92e16,6.18e16,6.53e15])

doPlot = True
doLoop = False

facA = 0.3239
#rptop_emp = np.sqrt((np_0+1/3*rp_ave*1e-4*dndy)/np_0)*rp_ave
rptop_emp = np.sqrt((np_0+facA*rp_ave*1e-4*dndy)/np_0)*rp_ave
rpbot_emp = np.sqrt((np_0-facA*rp_ave*1e-4*dndy)/np_0)*rp_ave

yoff_emp = 1/2*(rptop_emp-rpbot_emp)/(lhalf)

gset = np.where(np_0==2e16)[0]
yset = np.where(np_0==1e16)[0]
bset = np.where(np_0==1e17)[0]



facA_arr = np.array([facA])
fitqt_arr = np.zeros(1)
fitqb_arr = np.zeros(1)
"""
if doLoop:
    facA_arr = np.linspace(0.24,0.4,50)
    fitqt_arr = np.zeros(len(facA_arr))
    fitqb_arr = np.zeros(len(facA_arr))
    
for i in range(len(facA_arr)):
    
    if doLoop:
        facA = facA_arr[i]
        rptop_emp = np.sqrt((np_0+facA*rp_ave*1e-4*dndy)/np_0)*rp_ave
        rpbot_emp = np.sqrt((np_0-facA*rp_ave*1e-4*dndy)/np_0)*rp_ave
        
    rt_fit = np.polyfit(rptop_emp, rp_top, 1)
    rb_fit = np.polyfit(rpbot_emp, rp_bot, 1)
    
    if doPlot:
        plt.title("Emperical Guess for R_top")
        plt.xlabel(r'$R_{top,emp} \ (\mu m)$')
        plt.ylabel(r'$R_{top} \ (\mu m)$')
        plt.scatter(rptop_emp[gset],rp_top[gset],c='g', label="2e16 cm-3")
        plt.scatter(rptop_emp[yset],rp_top[yset],c='y', label="1e16 cm-3")
        plt.scatter(rptop_emp[bset],rp_top[bset],c='b', label="1e17 cm-3")
        plt.plot([30,110],[30,110],label="1:1", ls='--')
        #plt.plot(rptop_emp, rptop_emp*rt_fit[0]+rt_fit[1], ls='--')
        plt.grid(); plt.legend(); plt.show()
        print("Fit slope = ",rt_fit[0])
    
    if doPlot:
        plt.title("Emperical Guess for R_bot")
        plt.xlabel(r'$R_{bot,emp} \ (\mu m)$')
        plt.ylabel(r'$R_{bot} \ (\mu m)$')
        plt.scatter(rpbot_emp[gset],rp_bot[gset],c='g', label="2e16 cm-3")
        plt.scatter(rpbot_emp[yset],rp_bot[yset],c='y', label="1e16 cm-3")
        plt.scatter(rpbot_emp[bset],rp_bot[bset],c='b', label="1e17 cm-3")
        plt.plot([30,110],[30,110],label="1:1", ls='--')
        plt.grid(); plt.legend(); plt.show()
        print("Fit slope = ",rb_fit[0])
    
    if doPlot:
        plt.title("Error plot for Rtop and Rbot")
        plt.xlabel("Simulation Case")
        plt.ylabel("Error")
        plt.plot(-(rp_top-rptop_emp)/rp_top,label='top')
        plt.plot(-(rp_bot-rpbot_emp)/rp_bot,label='bot')
        plt.grid(); plt.legend(); plt.show()

    fitqt_arr[i] = np.abs(1-rt_fit[0])
    fitqb_arr[i] = np.abs(1-rb_fit[0])

if doLoop:
    plt.plot(facA_arr,fitqt_arr)
    plt.plot(facA_arr,fitqb_arr)
    plt.show()
    print("facA --> ", facA_arr[np.argmin(np.abs(fitqt_arr-fitqb_arr))])



facA_arr = np.array([facA])
fitqo_arr = np.zeros(1)
if doLoop:
    facA_arr = np.linspace(0.31,0.33,50)
    fitqo_arr = np.zeros(len(facA_arr))

for i in range(len(facA_arr)):
    
    if doLoop:
        facA = facA_arr[i]
        rptop_emp = np.sqrt((np_0+facA*rp_ave*1e-4*dndy)/np_0)*rp_ave
        rpbot_emp = np.sqrt((np_0-facA*rp_ave*1e-4*dndy)/np_0)*rp_ave
        yoff_emp = 1/2*(rptop_emp-rpbot_emp)/(lhalf)

    yoff_fit = np.polyfit(yoff_emp, yoff, 1)

    if doPlot:
        plt.title("Wave Offset Growth across Simulations")
        plt.xlabel(r'$\partial n / \partial y \ (cm^{-4})$')
        plt.ylabel(r'$y_{off} \, / \, \mu m$')
        plt.scatter(dndy[gset],yoff[gset],c='g', label="2e16 cm-3")
        plt.scatter(dndy[yset],yoff[yset],c='y', label="1e16 cm-3")
        plt.scatter(dndy[bset],yoff[bset],c='b', label="1e17 cm-3")
        #plt.plot([30,110],[30,110],label="1:1", ls='--')
        plt.grid(); plt.legend(); plt.show()
        
        plt.title("Emperical Guess for y_off")
        plt.xlabel(r'$y_{off,emp} \, / \, \mu m$')
        plt.ylabel(r'$y_{off} \, / \, \mu m$')
        plt.scatter(yoff_emp[gset],yoff[gset],c='g', label="2e16 cm-3")
        plt.scatter(yoff_emp[yset],yoff[yset],c='y', label="1e16 cm-3")
        plt.scatter(yoff_emp[bset],yoff[bset],c='b', label="1e17 cm-3")
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
        
    fitqo_arr[i] = np.abs(1-yoff_fit[0])

if doLoop:
    plt.plot(facA_arr,fitqo_arr)
    plt.show()
    print("facA --> ", facA_arr[np.argmin(fitqo_arr)])
"""


x = np.logspace(np.log10(min(dndy)),np.log10(max(dndy)),50)
y = np.zeros(len(x))
"""#We done with this, only need 1 scalar
power = 3
dndys_fit = np.polyfit(dndy, dndys, power)

#dndy_alt=np.array([0,0,0,0,8e17,2e17,4e18,2e18,4e17])
#dndys_alt = np.array([0,0,0,0,8.64e17,1.47e17,1.21e19,6.54e18,7.96e17])

#dndys_fit = np.polyfit(dndy_alt, dndys_alt, power)

dndys_fit[0]= -4.7757e-37#     -3.0789e-38
dndys_fit[1]= 2.7677e-18#     6.1568e-19
dndys_fit[2]= -0.4195#     1.0266
dndys_fit[3]= 6.4565e16#     2.5089e16

dndys_fit[3] = 0

print("Fit with a power of",power,":")
print("(dn/dy)_sheath = ")
for i in range(power+1):
    y = y+ np.power(x,power-i)*dndys_fit[i]
    print("    (dn/dy)^",power-i, " x " , dndys_fit[i]," +")

plt.title("Sheath Density Gradient")
plt.xlabel(r'$\partial n / \partial y \ (cm^{-4})$')
plt.ylabel(r'$\partial n / \partial y_s \ (cm^{-4})$')
plt.errorbar(dndy,dndys,yerr=dndys_std,xerr=None,linewidth = 0,elinewidth = 1,ecolor='k',barsabove = True)
plt.scatter(dndy[gset],dndys[gset],c='g', label="2e16 cm-3")
plt.scatter(dndy[yset],dndys[yset],c='y', label="1e16 cm-3")
plt.scatter(dndy[bset],dndys[bset],c='b', label="1e17 cm-3")
#plt.plot([30,110],[30,110],label="1:1", ls='--')
plt.grid(); plt.legend(); plt.show()
"""
eps = 8.854e-12
e = 1.602e-19

#empB = 2.2
#delta = 0.15

empB = 2.482
delta = 0.134

facB_arr = np.array([empB])
fitB_arr = np.zeros(1)
facD_arr = np.array([delta])
fitD_arr = np.zeros(1)

"""
#With LOW 1e17
dndy=np.array([8e17,2e17,2.5e16,2.5e15,2e18,4e17])
dndys = np.array([1.76e18,4.33e17,6.05e16,6.07e15,4.96e18,7.96e17])
dndys_std = np.array([1.18e17,4.44e16,4.02e16,4.65e15,3.15e18,1.88e17])
ering_norm = np.array([1.33e17,3.52e16,4.18e15,4.44e14,3.69e17,5.98e16])
ering_norm_std = np.array([7.86e15,3.91e15,1.64e15,4.28e14,6.18e16,6.53e15])
"""
"""
#With NO 1e17
dndy=np.array([8e17,2e17,2.5e16,2.5e15,4e17])
dndys = np.array([1.76e18,4.33e17,6.05e16,6.07e15,7.96e17])
dndys_std = np.array([1.18e17,4.44e16,4.02e16,4.65e15,1.88e17])
ering_norm = np.array([1.33e17,3.52e16,4.18e15,4.44e14,5.98e16])
ering_norm_std = np.array([7.86e15,3.91e15,1.64e15,4.28e14,6.53e15])
"""
"""
if doLoop:
    facB_arr = np.linspace(2.1,2.8,100)
    fitB_arr = np.zeros(len(facB_arr))
    
for i in range(len(facB_arr)):
    if doLoop:
        empB = facB_arr[i]
        sheath_fit = np.polyfit(dndy*empB, dndys, 1)
        sheath_qual = np.abs(1-sheath_fit[0])
        fitB_arr[i]=sheath_qual

    if doPlot:
        plt.title("Sheath Density Gradient")
        plt.xlabel(r'$\partial n / \partial y \ (cm^{-4})$')
        plt.ylabel(r'$\partial n / \partial y_s \ (cm^{-4})$')
        #plt.loglog(x,y,ls='--',label="Third Order Fit")
        #plt.loglog(dndy,dndy**1.02,ls='--',label="Power 1.02")
        plt.loglog(dndy,empB*dndy,ls='--',label="Scalar "+str(empB))
        plt.errorbar(dndy,dndys,yerr=dndys_std,xerr=None,linewidth = 0,elinewidth = 1,ecolor='k',barsabove = True)
        plt.scatter(dndy[gset],dndys[gset],c='g', label="2e16 cm-3")
        plt.scatter(dndy[yset],dndys[yset],c='y', label="1e16 cm-3")
        plt.scatter(dndy[bset],dndys[bset],c='b', label="1e17 cm-3")
        #plt.plot([30,110],[30,110],label="1:1", ls='--')
        plt.grid(); plt.legend(); plt.show()
        sheath_fit = np.polyfit(dndy*empB, dndys, 1)
        sheath_qual = np.abs(1-sheath_fit[0])
        print("Fit Quality: ",sheath_qual)
    
if doLoop:
    plt.plot(facB_arr,fitB_arr)
    plt.show()
    print("facB --> ", facB_arr[np.argmin(fitB_arr)])
"""
"""
ering_norm_righthand = e/2/eps*delta*((dndys-dndy)*100**4)
plt.title("Ering using AVERAGED dndy_sheath")
plt.xlabel(r'$\langle E_{ring}/R_p^2\rangle \ (V/m^3)$')
plt.ylabel("Eqn. 1 "+r'$(V/m^3)$')
plt.errorbar(ering_norm,ering_norm_righthand,yerr=ering_norm_std,xerr=None,linewidth = 0,elinewidth = 1,ecolor='k',barsabove = True)
plt.scatter(ering_norm[gset],ering_norm_righthand[gset],c='g', label="2e16 cm-3")
plt.scatter(ering_norm[yset],ering_norm_righthand[yset],c='y', label="1e16 cm-3")
plt.scatter(ering_norm[bset],ering_norm_righthand[bset],c='b', label="1e17 cm-3")
plt.plot([0,max(ering_norm)],[0,max(ering_norm)],label="1:1", ls='--')
plt.grid(); plt.legend(); plt.show()
"""
"""
#empB = 2.786
if doLoop:
    facD_arr = np.linspace(0.07,0.2,100)
    fitD_arr = np.zeros(len(facD_arr))
    
for i in range(len(facD_arr)):
    if doLoop:
        delta = facD_arr[i]
        ering_norm_righthand = e/2/eps*delta*(((empB-1)*dndy)*100**4)
        ering_fit = np.polyfit(ering_norm, ering_norm_righthand, 1)
        ering_qual = np.abs(1-ering_fit[0])
        fitD_arr[i]=ering_qual
    
    if doPlot:
        ering_norm_righthand = e/2/eps*delta*(((empB-1)*dndy)*100**4)
        plt.title("Ering using EMPIRICAL dndy_sheath")
        plt.xlabel(r'$\langle E_{ring}/R_p^2\rangle \ (V/m^3)$')
        plt.ylabel("Eqn. 2 "+r'$(V/m^3)$')
        plt.errorbar(ering_norm,ering_norm_righthand,xerr=ering_norm_std,yerr=None,linewidth = 0,elinewidth = 1,ecolor='k',barsabove = True)
        plt.scatter(ering_norm[gset],ering_norm_righthand[gset],c='g', label="2e16 cm-3")
        plt.scatter(ering_norm[yset],ering_norm_righthand[yset],c='y', label="1e16 cm-3")
        plt.scatter(ering_norm[bset],ering_norm_righthand[bset],c='b', label="1e17 cm-3")
        plt.plot([0,max(ering_norm)],[0,max(ering_norm)],label="1:1", ls='--')
        plt.grid(); plt.legend(); plt.show()
        ering_fit = np.polyfit(ering_norm_righthand, ering_norm, 1)
        ering_qual = np.abs(1-ering_fit[0])
        print("Fit Quality: ",ering_qual)

if doLoop:
    plt.plot(facD_arr,fitD_arr)
    plt.show()
    print("facD --> ", facD_arr[np.argmin(fitD_arr)])

"""
x = np.linspace(min(dndy),max(dndy),50)
sh = x * empB
y = e/2/eps*delta*(((empB-1)*x)*100**4)
ering_norm_righthand = e/2/eps*delta*(((empB-1)*dndy)*100**4)

fig, (ax0,ax1) = plt.subplots(nrows=2,ncols=1,sharex=True)
fig.set_size_inches(5,6)

pl0 = ax0.semilogx(x,sh,label="Empirical Sheath Gradient",ls='--')
ax0.errorbar(dndy,dndys,yerr=dndys_std,xerr=None,linewidth = 0,elinewidth = 1,ecolor='k',barsabove = True)
ax0.scatter(dndy[gset],dndys[gset],c='g', label=r'$\mathrm{n_0 = 2e16 \ cm^{-3}}$')
ax0.scatter(dndy[yset],dndys[yset],c='y', label=r'$\mathrm{n_0 = 1e16 \ cm^{-3}}$')
ax0.scatter(dndy[bset],dndys[bset],c='b', label=r'$\mathrm{n_0 = 1e17 \ cm^{-3}}$')
ax0.set_ylabel(r'$\langle\partial n / \partial y\rangle_{sh} \ \mathrm{(cm^{-3})}$')
ax0.legend()

E_si_to_cgs = (1/3)*10**(-4)
L_si_to_cgs = 100
factor = E_si_to_cgs/L_si_to_cgs/1e9
ering_norm = ering_norm*factor
ering_norm_std=ering_norm_std*factor
y = y*factor

pl0 = ax1.semilogx(x,y,label="Empirical Sheath Field", ls='--')
ax1.errorbar(dndy,ering_norm,yerr=ering_norm_std,xerr=None,linewidth = 0,elinewidth = 1,ecolor='k',barsabove = True)
ax1.scatter(dndy[gset],ering_norm[gset],c='g', label=r'$\mathrm{n_0 = 2e16 \ cm^{-3}}$')
ax1.scatter(dndy[yset],ering_norm[yset],c='y', label=r'$\mathrm{n_0 = 1e16 \ cm^{-3}}$')
ax1.scatter(dndy[bset],ering_norm[bset],c='b', label=r'$\mathrm{n_0 = 1e17 \ cm^{-3}}$')
ax1.set_ylabel(r'$\langle E_{ring}/R_p^2\rangle \ \mathrm{(GstatV/cm^3)}$')
ax1.set_xlabel(r'$(\partial n / \partial y)_i \ \mathrm{(cm^{-3})}$')
ax1.legend()

fig.subplots_adjust(hspace=0)

plt.show()