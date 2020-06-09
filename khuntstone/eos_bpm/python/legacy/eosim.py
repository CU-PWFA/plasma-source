'''
Module to run a full EOS-BPM signal simulation taking as input: THz electric
field(s), time domain of the field(s), and the probe pulse. Additionally can 
include separate positions for the drive and witness beams locations relative
to the crystal. Module also contains utility functions for performing parameter
scans and other optimizations. 

Created on Tue Jul 23 15:56:10 2019
@author: keenan
'''

# Imports and constants
from cycler import cycler;
import matplotlib.pyplot as plt;
import numpy as np;
from scipy.constants import c;
from scipy.signal import find_peaks;
# Custom modules
import phase_retard as pr; # For computing phase retardation of the probe
import thz_prop as prop; # For propagating THz fields in the crystal

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

def peak2peak(E, t, t_range):
    '''
    Function to find the spacing between the drive and witness beams peaks from 
    their electric field.
    Parameters:
    -----------
    E       : array_like
              The electric field
    t       : array_like
              Time array corresponding to E (s)
    t_range : array_like
              Array of [min, max] values of the offset between the beam. This is 
              used because sometimes the electric field is double peaked around 
              the drive beam so just finding the two highest values in the array
              is not helpful
    '''
    cutoff = 0;
    # Drive is already centered on 0
    peak1 = np.argmax(E);
    while True:
        cutoff += 0.05;
        peaks   = find_peaks(E, height = cutoff * max(E))[0];


        tdiff   = abs(t[peak1] - t[peaks[-1]]);
        if (t_range[0] <= tdiff and t_range[1] >= tdiff) or len(peaks) == 2:
            break
    return tdiff, peak1, peaks[-1];

def beam_peaks(drive_sig, wit_sig):
    '''
    Function to find the peak drive and peak witness signal
    Parameters:
    -----------
    drive_sig : array_like
                Array of the signal around the drive beam location
    wit_sig   : array_lik
                Array of the signal around the witness beam location
    Returns:
    --------
    drive_peaks : array_like
                  Array containing location of the drive peak signal
    wit_peaks   : array_like
                  Array containing location of the witness peak signal
    '''
    cutoff = 0.0;
    while True:
        cutoff += 0.05;
        drive_peaks = find_peaks(drive_sig,\
                        height = cutoff * max(drive_sig))[0];
        if len(drive_peaks) == 1:
            break;
        if len(drive_peaks) == 0:
            cutoff -= 0.5;
            drive_peaks = find_peaks(drive_sig, \
                          height = cutoff * max(drive_sig))[0];
            break;
    cutoff = 0.0;
    while True:
        cutoff += 0.5;
        wit_peaks  = find_peaks(wit_sig, height = cutoff * max(wit_sig))[0]
        if len(wit_peaks) == 1:
            break;
        if len(drive_peaks) == 0:
            cutoff -= 0.5;
            wit_peaks = find_peaks(wit_sig, \
                          height = cutoff * max(wit_sig))[0];
            break;
    return drive_peaks[0], wit_peaks[0];
def sig_peaks(E, t, gamma, t_plot):
    '''
    Function to compute and print (if wanted), beam peak to peak and signal
    peak to peak

    Parameters:
    -----------
    E      : array like
             THz field array or tuple of drive, witness, and combined fields
    t      : array_like
             Time array corresponding to E
    gamma  : array like
             array or tuple of phase(s) returned by sim (see below)
    t_plot : array_like
             Time array corresponding to gamma

    Returns:
    --------
    in_tdiff : float
               The peak to peak separation of the input THz pulse
    t_pol    : float
               The peak to peak separation of the signal at cross pol. 
    t_bal    : float
               The peak to peak separation of the signal using bal. detectros
    '''
    if type(E) == tuple:
        
        E_plot   = np.flip(E[2], axis = 0);
        te_plot  = (t - t[np.argmax(E_plot)]);
        t_wit    = te_plot[np.argmax(np.flip(E[1], axis = 0))]
        t_drive  = te_plot[np.argmax(np.flip(E[0], axis = 0))]
        in_tdiff = abs(t_drive - t_wit) 
        in_zdiff = in_tdiff * c * 1e6;
        print("Input peak to peak:", np.round(in_zdiff,2), "microns");

        pol_sig = np.sin(gamma[2] / 2)**2;
        pol_sig = pol_sig / max(pol_sig);

        bal_sig = np.sin(gamma[2]);
        bal_sig = bal_sig / max(bal_sig);

        # Get highest peak around drive / witness signals
        drive     = np.flip(E[0], axis = 0);
        wit       = np.flip(E[1], axis = 0);
        drive_ind = np.argmin(abs(t_plot - t_drive * 1e15))
        wit_ind   = np.argmin(abs(t_plot - t_wit * 1e15))

        drive_pol_sig = pol_sig[drive_ind-100:drive_ind+100]
        wit_pol_sig   = pol_sig[wit_ind-100:wit_ind+100]

        drive_peak, wit_peak = beam_peaks(drive_pol_sig, wit_pol_sig)
        drive_peak = drive_peak + drive_ind - 100;
        wit_peak   = wit_peak + wit_ind - 100;
        t_pol = abs(t_plot[drive_peak] - t_plot[wit_peak])
        z_pol = t_pol * 1e-15 * c;
        z_pol = np.round(z_pol * 1e6, 2);
        print("Crossed pol. peak to peak:", z_pol, "microns");

        drive_bal_sig = bal_sig[drive_ind-100:drive_ind+100]
        wit_bal_sig   = bal_sig[wit_ind-100:wit_ind+100]

        drive_peak, wit_peak = beam_peaks(drive_bal_sig, wit_bal_sig)
        drive_peak = drive_peak + drive_ind - 100;
        wit_peak   = wit_peak + wit_ind - 100;
        t_bal = abs(t_plot[drive_peak] - t_plot[wit_peak])
        z_bal = t_bal * 1e-15 * c;
        z_bal = np.round(z_bal * 1e6, 2);
        print("Balanced detectors peak to peak:", z_bal, "microns");

    return in_tdiff, t_pol, t_bal;
def sim(E, t, tau, ctype, d, probe, x = [0, 0], y = [1, 1], \
        plot = False, normed = False, verbose = False, plot_input = False, \
        plot_all = False):
    '''
    Simulates EOS-BPM signal
    Parameters:
    -----------
    E      : array like
             THz field array or tuple of drive, witness, and combined fields
    t      : array_like
             Time array corresponding to E
    tau    : array_like
             Time array corresponding to the probe
    ctype  : str
             The EO crystal (GaP or ZnTe)
    d      : float
             Crystal thickness (m)
    probe  : dictionary
             Dictionary containing keys 'y0', 'width', 'a_laser' with 
             corresponding values for central wavelength , wavelength range
             and angle of incidence (rad)
    x      : array_like, optional
             Array of horizontal offsets from crystal [drive, witness] (m)
    y      : array_like, optional
             Array of vertical offsets from crystal [drive, witness] (m)
    plot   : bool, optional
             Whether or not to plot the signal
    normed : bool, optional
             Whether or not to normalize the signal in plots
    '''
    # First compute the raw FFT of the fields
    FEr, f = prop.raw_field(E, t, plot = plot_all);
    # Propagate the fields through the given crystal
    Ec, tt, FEc, sf, t_shift = prop.cry_field(t, FEr, f, probe['y0'], ctype, d,\
        verbose = verbose, plot = plot_all);
    probe['t_shift'] = t_shift;
    # Compute the phase retardation
    gamma, t_plot = pr.phase_retard(Ec, tt, tau, probe, ctype, d, x = x, y = y);
    
    if verbose:
        dummy = sig_peaks(E, t, gamma, t_plot);

    
# TODO: The above print statements for a single THz field input
    if plot:
        if type(E) == tuple:
            # Just phase
            fig1, ax1 = makefig(x = 8, y = 6, xlab = 'time [fs]', \
                                ylab = r'$\Gamma$', \
                                title = 'Phase retardation', fs = 16, ts = 24);
            drive_sig = gamma[0]; 
            wit_sig   = gamma[1];
            sig       = gamma[2];
            ax1.plot(t_plot, sig, label = 'Combined');
            ax1.plot(t_plot, drive_sig, '--', label = 'Drive');
            ax1.plot(t_plot, wit_sig, '-.', label = 'Witness');
            ax1.legend();
            # Balanced detectors
            fig2, ax2 = makefig(x = 8, y = 6, xlab = 'time [fs]', \
                                ylab = 'Signal [AU]', \
                                title = 'Balanced detectors' \
                                 + r'$\propto \mathrm{sin}(\Gamma$)', 
                                 fs = 16, ts = 24);
            drive_sig = np.sin(gamma[0]); 
            wit_sig   = np.sin(gamma[1]);
            sig       = np.sin(gamma[2]);
            if normed:
                drive_sig = drive_sig / max(sig);
                wit_sig   = wit_sig / max(sig);
                sig       = sig / max(sig);
            ax2.plot(t_plot, sig, label = 'Combined');
            ax2.plot(t_plot, drive_sig,'--', label = 'Drive');
            ax2.plot(t_plot, wit_sig, '-.', label = 'Witness');
            ax2.legend();
            # Crossed polarizer
            fig3, ax3 = makefig(x = 8, y = 6, xlab = 'time [fs]', \
                                ylab = 'Signal [AU]', \
                                title = 'Cross pol.' \
                                 + r'$\propto \mathrm{sin}^2(\Gamma / 2)$', 
                                 fs = 16, ts = 24);
            drive_sig = np.sin(gamma[0] / 2)**2; 
            wit_sig   = np.sin(gamma[1] / 2)**2;
            sig       = np.sin(gamma[2] / 2)**2;
            if normed:
                drive_sig = drive_sig / max(sig);
                wit_sig   = wit_sig / max(sig);
                sig       = sig / max(sig);
            ax3.plot(t_plot, sig, label = 'Combined');
            ax3.plot(t_plot, drive_sig,'--', label = 'Drive');
            ax3.plot(t_plot, wit_sig, '-.', label = 'Witness');
            ax3.legend();

            # plot the input field:
            if plot_input:
                ax4 = ax1.twinx();
                ax5 = ax2.twinx();
                ax6 = ax3.twinx();
    
                ax4.tick_params(axis = 'y', labelcolor = 'red');
                ax5.tick_params(axis = 'y', labelcolor = 'red');
                ax6.tick_params(axis = 'y', labelcolor = 'red');
                
                E_plot = np.flip(E[2], axis = 0);
                te_plot = (t - t[np.argmax(E_plot)]) * 1e15
                ax4.plot(te_plot, E_plot / 1e6, '-r');
                ax4.set_ylabel(r'$E_{THz}$ [MV/m]', color = 'red');
    
                ax5.plot(te_plot, E_plot / 1e6, '-r');
                ax5.set_ylabel(r'$E_{THz}$ [MV/m]', color = 'red')
    
                ax6.plot(te_plot, E_plot / 1e6, '-r');
                ax6.set_ylabel(r'$E_{THz}$ [MV/m]', color = 'red')

                # Print the peak2peak separation of the input, and of the 
                # signal
            min_lim = max([min(t_plot), min(t * 1e15)]);
            max_lim = min([max(t_plot), max(t * 1e15)]);
            
            ax1.set_xlim([min_lim, max_lim]);
            ax2.set_xlim([min_lim, max_lim]);
            ax3.set_xlim([min_lim, max_lim]);
            ax3.set_xlim([min(t_plot), max(t_plot)]);
            plt.show();
    

        elif type(E) == np.ndarray or type(E) == list:
            fig1, ax1 = makefig(x = 8, y = 6, xlab = 'time [fs]', \
                                ylab = r'$\Gamma$', \
                                title = 'Phase retardation', fs = 16, ts = 24);
            ax1.plot(t_plot, gamma);

            fig2, ax2 = makefig(x = 8, y = 6, xlab = 'time [fs]', \
                                ylab = 'Signal [AU]', \
                                title = 'EOS-BPM Signal', fs = 16, ts = 24);
            cross = np.sin(gamma / 2)**2;
            bal   = np.sin(gamma);
            if normed:
                cross = cross / max(cross);
                bal   = bal / max(bal);
            ax2.plot(t_plot, bal, label = 'Balaced detectors');
            ax2.plot(t_plot, cross, '--',  label = 'Crossed pol.');
            ax2.legend();
            # Plot input field
            if plot_input:
                ax3 = ax1.twinx();
                ax3.tick_params(axis = 'y', labelcolor = 'red');
                ax3.plot(t * 1e15, E / 1e6, '-r');
                ax3.set_ylabel(r'$E_{THz}$ [MV / m]', color = 'red');
                
                ax4 = ax2.twinx();
                ax4.tick_params(axis = 'y', labelcolor = 'red');
                ax4.plot(t * 1e15, E / 1e6, '-r');
                ax4.set_ylabel(r'$E_{THz}$ [MV / m]', color = 'red');

            # Get plot range
            plot_min = max([min(t_plot), min(t * 1e15)]);
            plot_max = min([max(t_plot), max(t * 1e15)]);
            ax1.set_xlim([plot_min, plot_max]);
            ax2.set_xlim([plot_min, plot_max]);
            plt.show();          
    return gamma, t_plot;


