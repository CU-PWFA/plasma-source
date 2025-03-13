#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 11:46:48 2023

Loads a 3D density profile from the output of FourierRefraction and analyzes
the 2D planes and variance cuts of the profile.  Can also perform tanh and
Gaussian fits the the yz plane and x axis, respectively.

This copy of refractionanalysis is for my thesis, to plot panels across
multiple different simulations.  Just two cases here for two plots, but the rest
can be uncommented out if needed.

@author: chris
"""

import sys
sys.path.insert(0, "../")
from modules import ThreeDimensionAnalysis as ThrDim
from modules import TPLFocalLength as Foc
import numpy as np
import matplotlib.pyplot as plt

#'cuts' will plot density along axis and small variations
#'max_corrector' will shift the beam axis to the maximum density
#  poor for most cases-only interesting at high densities ~10e19
cuts=1
max_corrector=0
enforce_xoff=0#-60#60#-46#-31#for p2g8
#-11 for the spherical casea highp, -30 for highpz zoomed versions
calc_focal = 1
center_jet = 1 #set to 1 to not max correct on the jet axis

#'getfit' will fit the y and z axes to tanh profiles
#'fityx' uses getfit to add a Gaussian fit to the x axis
#'infinite_approx' compares getfit with simply assuming an infinite slab
focalfit = 1
getfit = 1
slanted = 1
fityx = 0
infinite_approx = 0

#size of window in micrometers
resize=1

caseSet=0
if caseSet == 0: #No Jet
    folder_arr = np.array([
            '/home/chris/Desktop/FourierPlots/cylindrical_ex_thesis_1e15/',
            '/home/chris/Desktop/FourierPlots/cylindrical_ex_thesis_3e16/',
            '/home/chris/Desktop/FourierPlots/cylindrical_ex_thesis_1e17/'
            ])
    n0_arr = np.array([1e15,3e16,1e17])
    
elif caseSet == 1: #Jet
    folder_arr = np.array([
            '/home/chris/Desktop/FourierPlots/cylindrical_ex_thesis_jet_3e16/',
            '/home/chris/Desktop/FourierPlots/cylindrical_ex_thesis_jet_1e17/',
            '/home/chris/Desktop/FourierPlots/cylindrical_ex_thesis_jet_1e18/'
            ])
    n0_arr = np.array([3e16,1e17,1e18])

# cylindrical_ex_thesis_xxx  and  cylindrical_ex_thesis_jet_xxx
# xxx = 1e15, 3e16, 1e17, and 1e18 for jet
#folder = '/home/chris/Desktop/FourierPlots/cylindrical_ex_thesis_jet_1e17/'
#n0_lens = 1e17
y_window = 400
z_window = 400
directory = 'case_1/'

plotSelect = 1 #0 for Transverse, 1 for Longitudinal

fig, axs = plt.subplots(3,1,sharex=True,figsize=(5,7))
cmap = plt.get_cmap('jet_r')
for i in range(len(folder_arr)):
    folder = folder_arr[i]
    n0_lens = n0_arr[i]

    path = folder+directory
    
    #Load densities and parameters, and morph Robert's output [beam, jet, laser]
    # into the size and orientation I use [laser, beam, jet]
    nplot = np.load(path+'finalDensity.npy')
    params = np.load(path+'params.npy').item()
    
    X = params['X']; Nx = params['Nx']
    Y = params['Y']; Ny = params['Ny']
    Z = params['Z']; Nz = params['Nz']
    
    den = ThrDim.RobertRoll(nplot)
    y = np.linspace(-X/2, X/2, Nx, False)
    z = np.linspace(-Y/2, Y/2, Ny, False)
    x = np.linspace(-Z/2, Z/2, Nz, False)
    
    if resize == 1:
        ywidth=int(round(len(y)/(y[0]/(-.5*y_window))))
        yi=int(round((len(y)-ywidth)/2))
        yf=yi+ywidth
        
        zwidth=int(round(len(z)/(z[0]/(-.5*z_window))))
        zi=int(round((len(z)-zwidth)/2))
        zf=zi+zwidth
        
        den=den[:,yi:yf,zi:zf]
        y=y[yi:yf]
        z=z[zi:zf]
    
    #No offset, unless we run max_corrector
    x_off=0; y_off=0; z_off=0
    """
    if max_corrector == 1:
        maximum = ThrDim.GetMaximumOffset(den)
        x_off=maximum[0]; y_off=maximum[1]; z_off=maximum[2]
        if (enforce_xoff != 0):
            x_off = + enforce_xoff
        if center_jet == 1:
            z_off = 0
        print('-Corrected beam position to maximum of '+str(maximum[3])+' e17cm^-3')
        print('-Corrected x plane at '+str(x[x_off + round(len(den[:,0,0])/2)])+' microns')
        print('-Corrected y plane at '+str(y[y_off + round(len(den[0,:,0])/2)])+' microns')
        print('-Corrected z plane at '+str(z[z_off + round(len(den[0,0,:])/2)])+' microns')
    """
    if cuts == 1:
        x_step=x[1]-x[0]
        y_step=y[1]-y[0]
        z_step=z[1]-z[0]
        both = True
        ylabel = 'Plasma Density ' r'$(\mathrm{cm^{-3}})$'
        if plotSelect == 0:
            den_plane_zx=np.transpose(den[:,round(len(y)/2),:])
            label=['Density along gas jet (Vary Laser Distance)',
                   'Radius from axis (microns)',
                   'ni (e17 cm^-3)',
                   'Offset in +/- x(microns)']
            #ThrDim.VarianceCut(den_plane_zx,z,x_off,3,1,x_step,label,True)
            
            data = den_plane_zx
            axis = z
            offset = x_off
            number = 3
            spacing = 1
            unit = x_step
            loc = 0
            
        
        if plotSelect == 1:
            den_plane_yz=den[round(len(x)/2)+x_off,:,:]
            #label=['Density along beam axis (Vary Jet Distance)',
            #       'Radius from axis (microns)',
            #       'ni (e17 cm^-3)',
            #       'Offset in z(microns)']
            #ThrDim.VarianceCut(den_plane_yz,y,z_off,5,2,z_step,label)
            #ThrDim.VarianceCut(den_plane_yz,y,z_off,5,-2,z_step,label)
            
            label=['Density along beam axis (Vary Jet Distance)',
                   r'$z\mathrm{\ (\mu m)}$',
                   r'$n_p \mathrm{\ (10^{17} cm^{-3}})$',
                   r'$\Delta y\mathrm{\ (\mu m)}$']
            #ThrDim.VarianceCut(den_plane_yz,y,z_off,5,5,z_step,label,True)
            #ThrDim.VarianceCut_Prod(den_plane_yz,y,z_off,5,2,z_step,label,True)
            #sys.exit();
            
            data = den_plane_yz
            axis = y
            offset = z_off
            number = 5
            spacing = 5
            unit = z_step
            loc = 2
    
        if plotSelect == 2:
            den_plane_yx=np.transpose(den[:,:,round(len(z)/2)+z_off])
            """#For the plus and minus individually
            label=['Density along beam axis (Vary Laser Distance)',
                   'Radius from axis (microns)',
                   'ni (e17 cm^-3)',
                   'Offset in x(microns)']
            ThrDim.VarianceCut(den_plane_yx,y,x_off,5,2,x_step,label)
            ThrDim.VarianceCut(den_plane_yx,y,x_off,5,-2,x_step,label)
            """
            
            label=['Density along beam axis (Vary Laser Distance)',
                   'Radius from axis (microns)',
                   'ni (e17 cm^-3)',
                   'Offset in +/- x(microns)']
            #ThrDim.VarianceCut(den_plane_yx,y,x_off,4,1,x_step,label,True)

            data = den_plane_yx
            axis = y
            offset = x_off
            number = 4
            spacing = 1
            unit = x_step
            loc = 2

        center = round(len(data[0,:])/2)+offset
        for j in np.arange(number):
            data_cut = data[:,center+j*spacing]
            col = cmap(float(j)/number)
            axs[i].plot(axis,data_cut,c=col,label=r'$\pm \ $'+"{:1.1f}".format(j*spacing*unit) + r'$\mathrm{\ \mu m}$')  
            if both:
                data_cut = data[:,center-j*spacing]
                axs[i].plot(axis,data_cut,c=col)
        axs[i].set_ylabel(ylabel)
        axs[i].legend(title=r'$n_0=$'+str(n0_arr[i])+r'$\mathrm{\ cm^{-3}}$',loc=loc)
        axs[i].grid()

if plotSelect == 0:
    axs[2].set_xlabel('Transverse Axis '+r'$(\mathrm{\mu m})$')
elif plotSelect == 1:
    axs[2].set_xlabel('Electron Beam Axis '+r'$(\mathrm{\mu m})$')
elif plotSelect == 2:
    axs[2].set_xlabel('Electron Beam Axis '+r'$(\mathrm{\mu m})$')
fig.tight_layout()
plt.show()
"""
sys.exit()

if focalfit == 1:
    den_vs_y=den[round(len(x)/2)+x_off,:,round(len(z)/2)]
    
    plt.plot(y, den_vs_y)
    plt.title("Density along beam axis")
    plt.xlabel(r'$\mathrm{Beam \ Axis \ [\mu m]}$')
    plt.ylabel(r'$\mathrm{Density \ [10^{17}cm^{-3}]}$')
    plt.grid(); plt.show()
    
    focal = Foc.Calc_Focus(den_vs_y,y)
    thick0 = Foc.Calc_Square_Lens(n0_lens, focal, Foc.gam_def)
    print(thick0,"um equivalent thickness")
    
    for plus in range(15):
        den_vs_y = den[round(len(x)/2)+x_off,:,round(len(z)/2) + plus]
        focal = Foc.Calc_Focus(den_vs_y,y)
        thick = Foc.Calc_Square_Lens(n0_lens, focal, Foc.gam_def)
        den_vs_ym = den[round(len(x)/2)+x_off,:,round(len(z)/2) - plus]
        focalm = Foc.Calc_Focus(den_vs_ym,y)
        thickm = Foc.Calc_Square_Lens(n0_lens, focalm, Foc.gam_def)
        print(thick,"um equivalent thickness at z= ", plus*z_step, " percent: ", (100*(thick-thick0)/thick0),"%")
        print(thickm,"um equivalent thickness at z= ", -plus*z_step, " percent: ", (100*(thickm-thick0)/thick0),"%")
    print();print("Now for horizontal variation:");print()
    for plus in (np.arange(15)*2):
        den_vs_y = den[round(len(x)/2)+x_off + int(plus),:,round(len(z)/2)]
        focal = Foc.Calc_Focus(den_vs_y,y)
        thick = Foc.Calc_Square_Lens(n0_lens, focal, Foc.gam_def)
        den_vs_ym = den[round(len(x)/2)+x_off - int(plus),:,round(len(z)/2)]
        focalm = Foc.Calc_Focus(den_vs_ym,y)
        thickm = Foc.Calc_Square_Lens(n0_lens, focalm, Foc.gam_def)
        print(thick,"um equivalent thickness at x= ", int(plus)*x_step, " percent: ", (100*(thick-thick0)/thick0),"%")
        print(thickm,"um equivalent thickness at x= ", -int(plus)*x_step, " percent: ", (100*(thickm-thick0)/thick0),"%")
    
#Fit the data to tanh in y, an elliptical tanh in yz, and Gaussian in x
if getfit == 1:
    horiz_offset = -6#-50
    den_vs_y=den[round(len(x)/2)+x_off+horiz_offset,:,round(len(z)/2)]
    
    if calc_focal == 1:
        Foc.Calc_Focus(den_vs_y,y)
        #den_ex = np.full((50),0.1)
        #y_ex = np.linspace(0,25,50)
        #Foc.Calc_Focus(den_ex,y_ex)
    
    #Calculate our guess by roughly estimating the length of the narrow waist
    dd = list(den_vs_y)
    Ltop = np.abs(y[dd.index(next(d for d in dd if d > 0))])
    print(Ltop)
    guess = [Ltop, Ltop/6., max(den_vs_y)]
    
    #Get the tanh paramters for the y axis by fitting den_vs_y
    fity = ThrDim.FitDataDoubleTanh(den_vs_y,y,guess,"Density vs y")
    
    #Get the tanh parameters for the z axis by fitting den_vs_z
    den_vs_z=den[round(len(x)/2)+x_off,round(len(y)/2),:]
    fitz = ThrDim.FitDataDoubleTanh(den_vs_z,z,[guess[0]/10.0,guess[1]*10.,guess[2]],"Density vs z")
    #plt.plot(z,den_vs_z)
    #fitz = [31.5,2,1]
    #plt.plot(z,ThrDim.DoubleTanh(fitz,z))
    #plt.show()
    #sys.exit()
    
    if slanted == 1:
        fitz_sl = ThrDim.FitDataSomething(den_vs_z,z,ThrDim.DoubleTanhSlant,[fitz[0],fitz[1],fitz[2],-0.01])
    
        comparePaper = 0
        if comparePaper == 1:
            fitz_old = [82.1345751991, 14.7258249296, 0.299171976932, -4.52975946e-04]
            plt.plot(z, den_vs_z)
            plt.plot(z, ThrDim.DoubleTanhSlant(fitz_old, z))
            plt.show()
    
    #Assume an elliptical tanh, plot the simulated and approximate yz planes
    # and take variance cuts through the 2D plane of their differences
    den_plane_yz=den[round(len(x)/2)+x_off,:,:]
    #p1 = ThrDim.Fit2DimTanh(den_plane_yz,y,z,[fity[0],fity[1],fity[2],fity[0]/fitz[0]])
    slant = 0
    if slanted == 1:
        slant = fitz_sl[3]
        print("slant: ",slant)
    difyz = ThrDim.Plot2DimDataTanh(den_plane_yz,y,z,fity,fitz,
                                    '(microns)','Plasma Density','e17(cm^-3)',slant)
    ThrDim.VarianceCut(np.transpose(difyz),y,0,5,5,z[1]-z[0],
                ['Plasma density difference along beam','Distance from axis (microns)',
                 'Density Difference e17(cm^-3)','Offset in z(microns)'],True)
        
    #Fit the x axis to a Gaussian and investigate the differences if we
    # multiply the tanh of the y axis to the Gaussian of the x axis
    if fityx == 1:
        den_vs_x=den[:,round(len(y)/2),round(len(z)/2)]
        fitx = ThrDim.FitDataGaussian(den_vs_x,x,[guess[2],guess[1]*40.,0],"Density vs x")
    
        den_plane_yx=np.transpose(den[:,:,round(len(z)/2)])
        difyx = ThrDim.Plot2DimTanhxGaussian(den_plane_yx,y,x,fity,fitx,
                    '(microns)','Plasma Density','e17(cm^-3)')
        ThrDim.VarianceCut(np.transpose(difyx),y,0,4,1,x[1]-x[0],
                ['Plasma density difference along beam','Distance from axis (microns)',
                 'Density Difference e17(cm^-3)','Offset in +/- x(microns)'],True)
    
    #Compare the yz plane simulated to a yz plane of an infinite slab
    if infinite_approx == 1:
        ThrDim.Plot2DimDataTanh(den_plane_yz,y,z,[fity[0],fity[1],fity[2],0],fitz,
                               '(microns)','Plasma Density','e17(cm^-3)')
"""