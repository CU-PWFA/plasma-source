'''
Module for computing the phase retardation undergone by a probe pulse scanning
a THz field in an EO crystal. 
'''

# site-packages imports
from cycler import cycler; # For making plots pretty
import matplotlib.pyplot as plt; # For plotting
import numpy as np;
from scipy.constants import c # Fundamental constants
import scipy.fftpack as sfft; # For FFTing fields
import scipy.interpolate as interpolate; # For interpolating finer eos signal

# custom modules
import eocry as ec;
import eosim as sim;
# colors for plotting

plot_colors = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', \
               '#DDCC77', '#CC6677', '#882255', '#AA4499'];
cy = cycler('color', plot_colors);

def makefig(x = 4, y = 3, xlab = '', ylab = '', title = '', fs = 12, ts = 12):
    '''
    Creates figure and axes handles.
    Parameters:
    -----------
    x     : int, optional
            Width of the figure (in), default = 4
    y     : int, optional
            Height of the figure (in), default = 3
    xlab  : str, optional
            x axis label, default = ''
    ylab  : str, optional
            y axis label, default = ''
    title : str, optional
            plot title, default = ''
    Returns:
    --------
    fig : object
          The figure handle
    ax  : object
          The axes handle
    '''
    fig = plt.figure(figsize = (x, y), dpi = 200);
    ax  = fig.gca();
    ax.set_prop_cycle(cy)
    ax.set_xlabel(xlab, fontsize = fs);
    ax.set_ylabel(ylab, fontsize = fs);
    ax.set_title(title, fontsize = ts, fontweight = 'bold');
    ax.tick_params(labelsize = 'large');
    return fig, ax;

def phase_retard(Ecs, tt, tau, probe, ctype, d, \
                 x = [0, 0], y = [1, 1], nslice = 100, plot = False):
    '''
    Function to compute the phase retardation experienced by the probing pulse
    in the EO crystal. 
    Parameters:
    -----------
    Ecs    : array_like
             Array of electric field in the crystal or tuple of arrays of
             drive, witness, and combined fields in the crystal
    tt     : array_like
             Time array corresponding to Ecs
    tau    : array_like
             Time array corresponding to probe pulse (s)
    probe  : dictionary
             Dictionary containing keys 'y0', 'width', 'a_laser' with 
             corresponding values for central wavelength , wavelength range
             and angle of incidence (rad)
    ctype  : str
             The type of EO crystal being used (GaP or ZnTe)
    d      : float
             The depth of the crystal
    x      : array_like, optional
             Array of horizontal offsets from crystal [drive, witness] (m)
    y      : array_like, optional
             Array of vertical offsets from crystal [drive, witness] (m)
    nslice : int
             Number of slices to divide the crystal into, default = 100
    plot   : bool
             Whether or not to plot the phase retardation and various signals
    '''

    y0      = probe['y0'];
    dy      = probe['width'];
    a_laser = probe['a_laser'];
    
    # Chirp the pulse
    y1 = y0 - (dy/2);
    y2 = y0 + (dy/2);
    # In Hz
    f1 = c / y1;
    f2 = c / y2;
    # Chirp parameter
    T = abs(tau[0] - tau[-1])
    k = (f2 - f1) / T;
    
    # Probe shift
    
    # Time array for calculating frequency must be positive
    t_freq = np.linspace(0, T, len(tau));
    f_arr  = f1 + k * t_freq;
    y_arr  = c / f_arr;
    
    # Create array of crystal slices
    j      = np.arange(1, nslice, 1);
    dz     = d / nslice;
    d_arr  = (j - 0.5) * dz;
    if type(Ecs) == tuple:

        Ec_drive  = Ecs[0];
        Ec_wit    = Ecs[1];
        Ec       = Ecs[2];
        
            
        # preallocate for loop
        gamma_drive = np.zeros(len(tau));
        gamma_wit   = np.zeros(len(tau));
        gamma       = np.zeros(len(tau));
        dz          = abs(d_arr[1] - d_arr[0]);
        
        n0               = ec.indref({'ctype' : ctype, 'lam' : y_arr});
        dummy, vg, dummy = ec.velocities({'ctype' : ctype, 'f' : f_arr});
        vg_eff           = vg * np.cos(probe['a_laser']);
        factor           = 2 * np.pi * n0**3 / y_arr;
        drive_factor     = (factor / 2) \
                            * np.sqrt(1 + 3 * x[0]**2 / (x[0]**2 + y[0]**2));
        wit_factor       = (factor / 2) \
                            * np.sqrt(1 + 3 * x[1]**2 / (x[1]**2 + y[1]**2));
        t_delay          = tau * 1e12;
        min_delay        = 0;
        max_delay        = 0;
        for j in range(len(d_arr)):
            t_interp = tt; 
            E_drive  = np.real(Ec_drive[:, j]);
            E_wit    = np.real(Ec_wit[:, j]);
            E        = np.real(Ec[:, j]);
            
            t_probe      = (d_arr[j] / vg_eff) * 1e12; # ps
            t_probe      = t_probe + probe['t_shift'] + t_delay;
            if max(t_probe) > max_delay:
                max_delay = max(t_probe);
            if min(t_probe) < min_delay:
                min_delay = min(t_probe);
            f_drive      = interpolate.interp1d(tt, E_drive, fill_value = 0.0, \
                                                bounds_error = False);
            f_wit        = interpolate.interp1d(tt, E_wit , fill_value = 0.0, \
                                                bounds_error = False);
            f            = interpolate.interp1d(tt, E, fill_value = 0.0, \
                                                bounds_error = False);
            gamma_drive += f_drive(t_probe) * dz;
            gamma_wit   += f_wit(t_probe) * dz;
            gamma       += f(t_probe) * dz;
        gamma_drive = gamma_drive * drive_factor;
        gamma_wit   = gamma_wit * wit_factor;
        gamma       = gamma * factor;
        
        
        # Recenter signals appropriately
        zero_ind = np.argmax(gamma);
        t_plot   = (tau - tau[zero_ind]) * 1e15; # fs
        if plot:
            fig1, ax1 = makefig(x = 8, y = 6, xlab = 't [fs]', \
                              ylab = r'$\Gamma$', title = 'Phase', \
                              fs = 16, ts = 24);
            ax1.plot(t_plot, gamma_drive, label = 'Drive');
            ax1.plot(t_plot, gamma_wit, label = 'Witness');
            ax1.plot(t_plot, gamma, '-.', label = 'Combined');
            ax1.legend();
            plt.show();
            
            fig2, ax2 = makefig(x = 8, y = 6, xlab = 't [fs]', \
                              ylab = r'Signal [AU]', \
                              title = \
                              r'Balanced detector' \
                              + r'$\propto \mathrm{sin}(\Gamma)$', \
                              fs = 16, ts = 24);
            ax2.plot(t_plot, np.sin(gamma_drive), label = 'Drive');
            ax2.plot(t_plot, np.sin(gamma_wit), label = 'Witness');
            ax2.plot(t_plot, np.sin(gamma), '-.', label = 'Combined');
            ax2.legend();
            plt.show();
            
            fig3, ax3 = makefig(x = 8, y = 6, xlab = 't [fs]', \
                              ylab = 'Signal [AU]', \
                              title = \
                              'Crossed pol.' \
                              + r'$\propto \mathrm{sin}^2(\Gamma / 2)$', \
                              fs = 16, ts = 24);
            ax3.plot(t_plot, \
                np.sin(gamma_drive / 2)**2 / max(np.sin(gamma / 2)**2),\
                       label = 'Drive');
            ax3.plot(t_plot, \
                np.sin(gamma_wit / 2)**2 / max(np.sin(gamma / 2)**2),\
                       label = 'Witness');
            ax3.plot(t_plot, np.sin(gamma / 2)**2 / max(np.sin(gamma / 2)**2),\
                     '-.', label = 'Combined');
            ax3.legend();
            plt.show();
        return (gamma_drive, gamma_wit, gamma), t_plot
    elif type(Ecs) == np.ndarray or type(Ecs) == list:
        Ec = Ecs;
        # preallocate for loop
        gamma       = np.zeros(len(tau));
        dz          = abs(d_arr[1] - d_arr[0]);
        n0               = ec.indref({'ctype' : ctype, 'lam' : y_arr});
        dummy, vg, dummy = ec.velocities({'ctype' : ctype, 'f' : f_arr});
        vg_eff           = vg * np.cos(probe['a_laser']);
        factor           = 2 * np.pi * n0**3 / y_arr;
        t_delay          = tau * 1e12;
        for j in range(len(d_arr)):
            E            = np.real(Ec[:, j]);
            t_probe      = (d_arr[j] / vg_eff) * 1e12; # ps
            t_probe      = t_probe + probe['t_shift'] + t_delay;
            t_interp     = tt;
            f            = interpolate.interp1d(t_interp, E, fill_value = 0, \
                                                bounds_error = False);
            #print(min(t_probe), max(t_probe), min(tt), max(tt))
            gamma       += f(t_probe) * dz;    
        gamma       = gamma * factor;
        # Recenter on 0
        zero_ind = np.argmax(gamma);
        t_plot = (tau - tau[zero_ind]) * 1e15; # fs
        if plot:
            fig, ax = makefig(x = 8, y = 6, xlab = 't [fs]', \
                              ylab = 'Signal [AU]');
            ax.plot(t_plot, gamma / max(gamma), label = r'$\Gamma$');
            ax.plot(t_plot, np.sin(gamma) / max(np.sin(gamma)),\
                    label = r'sin($\Gamma$)');
            ax.plot(t_plot, np.sin(gamma / 2)**2 / max(np.sin(gamma / 2)**2),\
                    label = r'sin$^2(\Gamma / 2)$');
            ax.legend();
            plt.show();
        return gamma, t_plot;
    else:
        print("ERROR: incorrect input electric field, returning null")
        return; 
