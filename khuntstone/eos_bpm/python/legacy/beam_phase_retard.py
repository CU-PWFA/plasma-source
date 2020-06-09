"""

Module to simulate detector signal from two THz pulses in an EO crystal. 
Created on Fri May 31 14:41:40 2019

@author: keenan
"""
# imports and constants
import copy;
from cycler import cycler
import matplotlib.pyplot as plt;
import numpy as np;
from scipy.constants import c;
import scipy.interpolate as interpolate;
import eocry as ec;
import two_beam as tb;

# Colors for plotting.
plot_colors = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', \
               '#DDCC77','#CC6677', '#882255',  '#AA4499'];
cy = cycler('color', plot_colors);


def phase_retard(params):
    
    
    '''
    Function to compute the phase retardation throughout a given EO crystal. 
    Parameters:
    -----------
    params : dictionary
             The params dictionary should contain
                ctype   : str
                          Specifies the EO crystal, either GaP or ZnTe
                t       : array_like
                          Time array for computing THz pulse in s.
                tau     : array_like
                          Delay time for the laser pulse in s. 
                d       : float
                          Crystal thickness in m.
                nslice  : int
                          The number of slices to divide the crystal into.
                beam1   : dictionary
                          A dictionary containing 'Q':charge (C), 
                          'x' : horizontal offset from crystal (m), 
                          'y' : vertical 'sig_t': temporal rms length (s)
                          for the primary bunch.
                beam2   : dictionary
                          A dictionary containing the same keys as beam1, with 
                          values for the trailing bunch.
                t_off   : float
                          The longitudinal temporal offset of the bunches. 
                probe   : dictionary, optional
                          A dictionary containing
                          'y0': float, the central wavelenth in m 
                          'z_match' : float, the depth at which the probe is 
                                      aligned with the peak of the drive pulse 
                          'a_laser' : float, the angle the laser makes with the 
                                      EO crystal in rads 
                          'chirp' : bool, whether or not we are using a chirped
                                    pulse for spectral encoding. 
                          'width' : float, the width of the chirped pulse in nm
                        
                plot    : bool, optional
                          Whether or not to plot the phase retardation.
                pockel  : bool, otpional
                          Whether or not to account for the Pockel effect when
                          propagating the THz pulse. default = True.
                ref     : bool, optional
                          Whether or not to include reflection of the THz pulse
                          default = False.   
    Returns:
    --------
    gamma_drive : array_like
                  The phase shift from the drive beam
    gamma_wit   : array_like
                  The phase shift from the witness beam
    gamma       : The combined phase shift
    
    NOTE: Since the drive, witness, and combined beams are propagated separately
    gamma != gamma_drive + gamma_wit unless the drive and witness beams are 
    *not* offset transversely.
    '''
    # Extract params
    ctype = params['ctype']
    probe = params['probe'];  
    tau   = params['tau'];
    det   = params['det'];
    y0    = probe['y0']
    beam1 = params['beam1'];
    beam2 = params['beam2'];
    x1    = beam1['x'];
    y1    = beam1['y'];
    x2    = beam2['x'];
    y2    = beam2['y'];
    #print("phase retard t_off", params['t_off']);
    if bool(params.get('plot')):
        plot = params['plot'];
    else:
        plot = False;                    
    if bool(probe.get('chirp')):
        chirp = probe['chirp'];
    else:
        chirp = False;

       
    # Get the temporal profiles of the pulses in the crystal at all slices
    cry_params = params; cry_params['plot'] = False; 
    Ec_drive, Ec_wit, Ec, tt, FEc_drive, FEc_wit, FEc, sf, d_arr = \
                                                  tb.thz_field_cry(cry_params);
       
    if chirp:
        # Get probe wavelength as a function of time
        del_y = probe['width'];
        # Down-chirp (frequency decreases in time)
        y1    = y0 - (del_y / 2);
        y2    = y0 + (del_y/2);
        # in Hz
        f1    = c / y1;
        f2    = c/y2;
        # chirp parameter
        T = abs(tau[0] - tau[-1])
        k = (f2 - f1) / T;
        
        # time array for calculating frequency must be positive
        # still corresponds to frequency as function of delay. 
        t_freq = np.linspace(0, T, len(tau));
        f_arr  = f1 + k * t_freq;
        y_arr  = c / f_arr;

        # preallocate for loop
        gamma_drive = np.zeros(len(tau));
        gamma_wit   = np.zeros(len(tau));
        gamma       = np.zeros(len(tau));
        factor      = np.zeros(len(tau));
        dz          = abs(d_arr[1] - d_arr[0]);
        
        n0               = ec.indref({'ctype' : ctype, 'lam' : y_arr});
        dummy, vg, dummy = ec.velocities({'ctype' : ctype, 'f' : f_arr});
        vg_eff           = vg * np.cos(probe['a_laser']);
        factor           = 2 * np.pi * n0**3 / y_arr;
        drive_factor     = (factor / 2) \
                                     * np.sqrt(1 + 3 * x1**2 / (x1**2 + y1**2));
        wit_factor       = (factor / 2) \
                                     * np.sqrt(1 + 3 * x2**2 / (x2**2 + y2**2));
        t_delay          = tau * 1e12;
        print("PROBE:", probe['t_shift'])
        for j in range(len(d_arr)):
            E_drive = np.real(Ec_drive[:, j]);
            E_wit   = np.real(Ec_wit[:, j]);
            E       = np.real(Ec[:, j]);
            
            t_probe      = (d_arr[j] / vg_eff) * 1e12; # ps
            t_probe      = t_probe + probe['t_shift'] + t_delay;
            f_drive      = interpolate.interp1d(tt, E_drive);
            f_wit        = interpolate.interp1d(tt, E_wit);
            f            = interpolate.interp1d(tt, E);
            gamma_drive += f_drive(t_probe) * dz;
            gamma_wit   += f_wit(t_probe) * dz;
            gamma       += f(t_probe) * dz;
                
        gamma_drive = gamma_drive * drive_factor;
        gamma_wit   = gamma_wit * wit_factor;
        gamma       = gamma * factor;
            
    else:
        # Get parameters of the single delta-like probe:
        f_opt = [c / y0]; # Hz;
        n0 = ec.indref({'ctype' : ctype, 'lam' : y0})
        dummy, vg, dummy = ec.velocities({'ctype' : ctype, 'f' : f_opt})
        
        vg_eff = vg * np.cos(probe['a_laser']);
        # preallocate
        gamma_drive = np.zeros(len(tau));
        gamma_wit   = np.zeros(len(tau));
        gamma       = np.zeros(len(tau));
        dz          = abs(d_arr[1] - d_arr[0]);
        factor = 2 * np.pi * n0**3 / y0;
        drive_factor     = (factor / 2) \
                                     * np.sqrt(1 + 3 * x1**2 / (x1**2 + y1**2));
        wit_factor       = (factor / 2) \
                                     * np.sqrt(1 + 3 * x2**2 / (x2**2 + y2**2));
        # tt is in ps, need to convert
        t_delay = tau * 1e12;
        for j in range(len(d_arr)):
            E_drive = np.real(Ec_drive[:, j]);
            E_wit   = np.real(Ec_wit[:, j]);
            E       = np.real(Ec[:, j]);
            t_probe      = (d_arr[j] / vg_eff) * 1e12; # ps
            t_probe      = t_probe + probe['t_shift'] + t_delay;
            f_drive      = interpolate.interp1d(tt, E_drive);
            f_wit        = interpolate.interp1d(tt, E_wit);
            f            = interpolate.interp1d(tt, E);
            
            gamma_drive += f_drive(t_probe) * dz;
            gamma_wit   += f_wit(t_probe) * dz;
            gamma       += f(t_probe) * dz;
        gamma_drive = gamma_drive * drive_factor;
        gamma_wit   = gamma_wit * wit_factor;
        gamma       = gamma * factor;
    if plot:
        fig = plt.figure(figsize = (6,6), dpi = 200);
        ax  = fig.gca();
        ax2 = ax.twinx()
        ax2.set_ylabel(r'$E_{thz}$', color = 'red');
        ax2.tick_params(axis = 'y', labelcolor = 'red');
        ax.set_prop_cycle(cy);
        # recenter on 0
        zero_ind = np.argmax(gamma);
        t_plot = (tau - tau[zero_ind]) * 1e15;
        if det == 'bal':
            combo_plot = np.sin(gamma);
            drive_plot = np.sin(gamma_drive);
            wit_plot   = np.sin(gamma_wit);
        elif det == 'cross':
            combo_plot = np.sin(gamma / 2)**2;
            drive_plot = np.sin(gamma_drive / 2)**2;
            wit_plot   = np.sin(gamma_wit / 2)**2;
        else:
            combo_plot = gamma;
            drive_plot = gamma_drive;
            wit_plot   = gamma_wit;
        ax.plot(t_plot, combo_plot, '-', label = 'Combined');
        ax.plot(t_plot, drive_plot, '--', label = 'Drive');
        ax.plot(t_plot, wit_plot, '-.', label = 'Witness');
        # Get the original pulse for plotting
        raw_params = {'beam1' : params['beam1'], 'beam2' : params['beam2'],
                      't_off' : - params['t_off'], 't' : t_plot * 1e-15}
        Er1, Er2, Er, t, FEr1, FEr2, FEr, f = tb.thz_field_raw(raw_params);
        ax2.plot(t_plot, (Er / max(Er)) * max(combo_plot), '-r');
        ax.legend()
        ax.set_xlabel(r'$\tau$ [fs]');
        if det == 'bal':
            ax.set_ylabel(r'sin$\Gamma$');
        elif det == 'cross':
            ax.set_ylabel(r'sin$^2\frac{\Gamma}{2}$');
        else:
            ax.set_ylabel(r'$\Gamma$ [rad]');
        if params['save']:
            plt.savefig(params['sname']);
        plt.show()                                          
    zero_ind = np.argmax(gamma)
    return gamma_drive, gamma_wit, gamma, (tau - tau[zero_ind]) * 1e15;                                                       
                                                       
                                                       
                                                       
                                            
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
