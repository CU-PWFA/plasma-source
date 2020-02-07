
"""
Created on Tue Feb 27 11:33:11 2018

@author: keenan
"""

import numpy as np
import sys
import KH_PlasmaProp as PProp
from matplotlib import pyplot as plt
sys.path.insert(0, "../litos")
import beam_ana as ba
import scipy.stats as stats
debug = 1
prop = 1
def PlotPropagation(ebeam, vbeam, plasma):
    nstep     = len(ebeam)
    s         = np.zeros(nstep)
    beta      = np.zeros(nstep)
    v_beta    = np.zeros(nstep)
    rms_x_eps = np.zeros(nstep)
    J_kurt    = np.zeros(nstep)
    frac      = 1.00
    
    for i in range(0,nstep):
        ebeam_rms = ba.calc_ebeam_rms(ebeam,i,frac)
        s[i]     = ebeam[i]["s"]
        beta[i] = ebeam_rms["x_beta"]/(1e-2)
        v_beta[i] = vbeam[i]["beta"]/(1e-2)
        rms_x_eps[i] = ebeam_rms["x_eps"]/(1e-6)
        
        [u,v] = ba.real2norm_coords(ebeam[i]["x"],ebeam[i]["xp"],\
                                ebeam_rms["x_beta"],ebeam_rms["x_alpha"])
        J = (u**2+v**2)/2
        J_kurt[i] = stats.kurtosis(J,0,False,True)
    
    figA, (ax1, ax3) = plt.subplots(2, sharex=True, sharey=False)
    
    ax1.plot(s,v_beta/10,color='b',linestyle='dashed')
    ax1.plot(s,beta,color='b',linestyle='solid')
    ax1.set_xlim([0.5,3.0])
    ax1.set_ylim([0,4.0])
    ax1.set_ylabel(r'$\beta$ [cm]',color='b')
    ax1.tick_params('y',colors='b')
    
    npl = plasma["npl"]/plasma["bulk"]["npl0"]
    
    ax2  = ax1.twinx()
    ax2.plot(s,npl,color='g',linestyle='solid')
    ax2.set_ylabel(r'$n_p/n_{p,0}$',color='g')
    ax2.tick_params('y',colors='g')
    ax2.set_ylim([0,1.4])
    ax2.text(0.50, 0.80, r'$n_{p,0} = %2.1e$'%plasma["bulk"]["npl0"],
            verticalalignment='center', horizontalalignment='center',
            transform=ax2.transAxes,
            color='green', fontsize=12)
    
    BB = PProp.CalcBmag(ebeam, plasma)
    print('Bmag: ', BB)
    ax3.plot(s,rms_x_eps/rms_x_eps[0],color='k',linestyle='-')
    ax3.plot(s,BB*np.ones(len(s)),color='k',linestyle='-.')
    ax3.set_ylabel(r'$\varepsilon_n/\varepsilon_{n,0}$',color='k')
    ax3.tick_params('y',colors='k')
    ax3.set_xlim([0.5,3.0])
    ax3.set_xlabel('z [m]')
    ax3.set_ylim([0.975,1.075])
    
    ax4 = ax3.twinx()
    ax4.plot(s,J_kurt/J_kurt[0],color='r',linestyle='--')
    ax4.set_ylabel(r'$K_J/K_{J,0}$',color='r')
    ax4.tick_params('y',colors='r')
    ax4.set_ylim([0.975,1.075])
    """
    #These act goofy
    xlabel_locs = [0.5,1.0,1.5,2.0,2.5,3.0]
    xlabels = [0,0.5,1.0,1.5,2.0,2.5]
    plt.xticks(xlabel_locs, xlabels)
    
    ax1.set_xlim([0.5,3.0])
    ax3.set_xlim([0.5,3.0])
    #end goofs
    """
    figA.tight_layout()
    figA.subplots_adjust(hspace=0)
    
    plt.show()


# Match waist location and half-width 
waist0 = -.272;
hw0    = .119;
ds    = 1e-2;
eps = .1
#num should be much higher, but 10 can give a good rough estimation
waist_arr = np.arange(waist0 - eps, waist0 + eps, ds, dtype = 'float');
hw_arr    = np.arange(hw0 - eps, hw0 + eps + ds, ds, dtype = 'float');



params  = PProp.ReturnDefaultParams()
params['npart'] = 0;
#params['L_up']  = 0;
params['L_dn']  = 0;

def bmagLoop(waist_arr, hw_arr, bmag_image):
    for i in range(len(waist_arr)):
        if debug == 1: print(i, 'of', len(waist_arr))
        for j in range(len(hw_arr)):
            waist   = waist_arr[i]
            hw      = hw_arr[j]
            params['hw_up'] = hw; params['hw_dn'] = hw; params['waist'] = waist;
            params['s_w'] = params['L_up'] + waist;
            twiss   = PProp.CallMakeTwiss(params)
            parts   = PProp.CallMakeParts(twiss, params)
            ebeam   = PProp.CallMakeBeam(twiss, parts, params)
            ebeam   = PProp.PropagateBackwards(ebeam, params)
            plasma  = PProp.MakeBulkPlasma(params)
            ebeam   = PProp.PropagatePlasma(ebeam, plasma)
            Bmag = PProp.CalcBmag(ebeam, plasma)
            #if debug == 1: print('Bmag: ',Bmag,'z: ',waist,'hw: ',hw)
            bmag_image[i][j] = Bmag
    return bmag_image
#print('initial Bmag: ',Bmag_init)
while True:
    print(ds)
    waist_arr = np.arange(waist0 - eps, waist0 + eps, ds, dtype = 'float');
    hw_arr    = np.arange(hw0 - eps, hw0 + eps + ds, ds, dtype = 'float');
    bmag_image = np.zeros((len(waist_arr),len(hw_arr)))
    bmag_image = bmagLoop(waist_arr, hw_arr, bmag_image)
    i_Bmin_x = np.argmin(np.min(bmag_image,1))
    i_Bmin_y = np.argmin(np.min(bmag_image,0))
    waist0 = waist_arr[i_Bmin_x]
    hw0    = hw_arr[i_Bmin_y]
    ds = ds/10;
    eps = eps/10
    if ds < 1e-6:
        print('Final coarser run')
        ds = 1e-3
        waist_arr = np.arange(waist0 - .1, waist0 + .1, ds, dtype = 'float')
        hw_arr    = np.arange(hw0 - .1, hw0 + .1, ds, dtype = 'float')
        bmag_image = np.zeros((len(waist_arr),len(hw_arr)))
        bmag_image = bmagLoop(waist_arr, hw_arr, bmag_image)
        break
minloc = PProp.PlotContour(bmag_image, waist_arr, hw_arr, \
                               r'$z_{\beta*}$ [m]', r'$\sigma_{hw}$ [m]')


mwaist = minloc[0]; mhw = minloc[1]

params = PProp.ReturnDefaultParams(hwup_change = mhw, waist_change = mwaist)
#params['L_up']  = 1.5;
if prop == 1:
    twiss           = PProp.CallMakeTwiss(params)
    parts           = PProp.CallMakeParts(twiss, params)
    ebeam0          = PProp.CallMakeBeam(twiss, parts, params)
    ebeam0          = PProp.PropagateBackwards(ebeam0, params)
    vbeam0          = ebeam0.copy();
    plasma0         = PProp.MakeBulkPlasma(params)
    ebeam0          = PProp.PropagatePlasma(ebeam0, plasma0)
    vbeam0          = PProp.PropagateVirtual(vbeam0, plasma0)
    PlotPropagation(ebeam0, vbeam0, plasma0)

