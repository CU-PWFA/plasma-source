'''
Module for propagating a THz pulse (or pulses) Through an EO crystal
'''


from cycler import cycler; # For making plots pretty
import matplotlib.pyplot as plt; # For plotting
import numpy as np;
from scipy.constants import c, epsilon_0 # Fundamental constants
import scipy.fftpack as sfft; # For FFTing fields


# Custom modules
import eocry as ec;

eps0 = epsilon_0;

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
# Propagation functions:
def prop_field(FT, A, r41, f_Hz, n, kappa, d, t, nslice):
    '''
    Propagates electric field through an EO crystal
    Parameters:
    -----------
    FT    : array_like
            The FFT of the electric field
    A     : array_like
            Transmission coefficient of the EO crystal as a function of 
            frequency
    r41   : array_like
            Pockel coefficient of the EO crysal as a function of frequency
    f_Hz  : array_like
            Frequency array (Hz)
    n     : array_like
            Index of refraction of the crystal as a function of frequency
    kappa : array_like
            Attenuation factor of the crystal as a function of frequency
    d     : float
            Crystal thickness in m
    t     : array_like
            Time array corresponding to the field 
    '''

    # FFT params
    L = len(t);

    # Create array of depths
    j      = np.arange(1, nslice, 1);
    dz     = d / nslice;
    d_arr  = (j - 0.5) * dz;
    # Preallocate for loop
    FEc             = np.zeros((len(f_Hz), len(d_arr)), \
                      dtype = 'complex128');
    Ec              = np.zeros((len(f_Hz), len(d_arr)), \
                      dtype = 'complex128');
    for i in range(len(d_arr)):
        FEc[:,i] = A * FT \
                     * r41 * np.exp((1j * 2 * np.pi * f_Hz / c) \
                     * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                     * d_arr[i]);
    
        Ec[:,i]  = sfft.ifft(FEc[:,i]) * L;
        #Ec[:, i] = np.flip(Ec[:, i], axis = 0); 
    return FEc, Ec;

def prop_ref(FT, A, r41, A_ref, f_Hz, n, kappa, d, t, nslice):
    '''
    Propagates electric field reflection through an EO crystal
    Parameters:
    -----------
    FT    : array_like
            The FFT of the electric field
    A_ref : array_like
            Reflection coefficient of the EO crystal as a function of frequency
    f_Hz  : array_like
            Frequency array (Hz)
    n     : array_like
            Index of refraction of the crystal as a function of frequency
    kappa : array_like
            Attenuation factor of the crystal as a function of frequency
    d     : float
            Crystal thickness in m
    t     : array_like
            Time array corresponding to the field 
    '''
    L = len(t);
    # Create array of depths
    j      = np.arange(1, nslice, 1);
    dz     = d / nslice;
    d_arr  = (j - 0.5) * dz;

    k         = 2 * np.pi * f_Hz * n / c; 
    alpha     = 2 * np.pi * f_Hz * kappa / c;
    ref_phase = np.exp(1j * 2 * k * d) * np.exp(-2 * alpha * d);

    # Preallocate for loop
    F2    = np.zeros((len(f_Hz), len(d_arr)), dtype = 'complex128');
    E_ref = np.zeros((len(f_Hz), len(d_arr)), dtype = 'complex128');
    for i in range(len(d_arr)):
        F2[:, i]       = A * FT * r41 * A_ref**2 * ref_phase  \
                           * np.exp((1j * 2 * np.pi * f_Hz / c) \
                           * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                           * d_arr[i]);                                              
        E_ref[:, i]       = sfft.ifft(F2[:, i]) * L;
        E_ref[:, i]       = np.flip(E_ref[:, i], axis = 0);
    return F2, E_ref;

# Verbose functions
def print_vel(ctype, f_Hz, verbose = False):
    # Calculates and prints velocities
    v_ph, v_g, dum = ec.velocities({'ctype' : ctype, 'f' : f_Hz});
    v_ph_c = np.round(np.nanmean(v_ph / c), 2);
    v_g_c  = np.round(np.nanmean(v_g / c), 2);

    if verbose:
        print(f'THz average phase velocity = {v_ph_c}c');
        print(f'THz average group velocity = {v_g_c}c');

def shift_probe(y0, Ec, tt, d, nslice, ctype, a_laser = 0, verbose = False):
    # Create array of depths
    j      = np.arange(1, nslice, 1);
    dz     = d / nslice;
    d_arr  = (j - 0.5) * dz;
    f_opt  = np.array([c / y0]);
    dum, vg_opt, dum = \
                      ec.velocities({'ctype' : ctype, 'f' : f_opt});
    vg_eff = vg_opt * np.cos(a_laser);
    # Add shift so probe is aligned to z_match
    #match_ind = np.argmin(abs(d_arr - 9.5e-6));
    match_ind = 0;
    t_peak  = tt[np.argmax(np.real(Ec[:, match_ind]))];
    t_shift = t_peak - ((d_arr[match_ind] / vg_eff) * 1e12);
    if verbose:
        print("Effective probe velocity =", \
              str(np.round(vg_eff / c, 2)), "c");
        print("Probe shift:" + str(np.round(t_shift, 2)) + "ps");
    return t_shift, vg_eff;

def raw_field(E, t, plot = False):
    '''
    Computes the FT of the input field and (as a check) its IFFT.

    Parameters:
    ----------
    E : array_like
        The input THz electric field or a 2d tuple containing arrays for 
        the drive THz field and the witness THz field
        t : array_like
        The time array corresponding to the fields (s)
    plot : bool, optional
           Whether or not to plot the input field, its FFT, and the IFFT
    Returns:
    -------
    FEr : array_like
          Array of the input field FFT or tuple of the FFTs for the drive
          and witness beam
    f   : array_like
          Frequency corresponding to FEr in THz
    '''
    def nextpow2(i):
        n = 1;
        while 2**n < i : n += 1;
        return n;
    def get_FFTs(E, t):
        Ts = abs(t[1] - t[0]) # sampling time interval (s)
        Fs = 1/Ts; # Sampling frequency (Hz)
        L  = len(t); # Number of samples

        def nextpow2(i):
            n = 1;
            while 2**n < i : n += 1
            return n
        # Get the FFT
        NFFT   = 2**(nextpow2(L) + 1);
        f      = Fs * np.linspace(0, 1, NFFT)*1e-12; # THz frequencies
        f_plot = (Fs/2) * np.linspace(0, 1, int(NFFT / 2 + 1)) * 1e-12; 
        FEr    = sfft.fft(E, NFFT) / L;
        
        # As a general check calculate the IFFT
        Nsamp = 2 * (len(f_plot) - 1);
        fmax  = f_plot[-1]*1e12; # Nyquist freqency Hz
        tt    = np.arange(0, Nsamp - 1) * (1 / (2 * fmax)); # s, time axis
        IFER  = sfft.ifft(FEr) * L;
        
        # Remove extra zeros
        tdiff = len(tt) - len(t);
        tt    = tt[0:-1 - tdiff + 1];
        IFER  = IFER[0:-1 - tdiff];
        
        # recenter time axis on '0', only works if input t is symmetric
        tt    = tt - tt[np.argmax(E)]
        tt    = tt * 1e15; # ps
        
        return FEr, f, IFER, tt, f_plot
    if type(E) == tuple:
        E_drive = E[0];
        E_wit   = E[1];
        E_both  = E_drive + E_wit;
        FE_drive, f, IFT_drive, tt, f_plot = get_FFTs(E_drive, t);
        FE_wit, f, IFT_wit, tt, f_plot = get_FFTs(E_wit, t);
        FE, f, IFT, tt, f_plot = get_FFTs(E_both, t);
        # Plot if requested
        if plot:
            NFFT = 2**(nextpow2(len(t)) + 1);
            fig1, ax1 = makefig(xlab = 't [fs]', \
                                ylab = 'E [V/m]', \
                        title = 'Input Field');
            ax1.plot(t * 1e15, E_drive, '-', label = 'Drive');
            ax1.plot(t * 1e15, E_wit, '-', label = 'Witness');
            ax1.plot(t * 1e15, E_both, '-.', label = 'Combined');
            ax1.legend();

            fig2, ax2 = makefig(xlab = 'f [THz]', \
                                ylab = r'$F[E_r]$', \
                        title = 'FFT');
            FE_drive_plot = 2 * abs(FE_drive[1:int(NFFT / 2 + 2)]);
            FE_wit_plot   = 2 * abs(FE_wit[1:int(NFFT / 2 + 2)]);
            FE_plot       = 2 * abs(FE[1:int(NFFT / 2 + 2)]);
            ax2.plot(f_plot, FE_drive_plot,'-', label = 'Drive');
            ax2.plot(f_plot, FE_wit_plot, '-', label = 'Witness');
            ax2.plot(f_plot, FE_plot, '-.', label = 'Combined'); 
            ax2.set_xlim([0, 20]);
            ax2.legend();

            fig3, ax3 = makefig(xlab = 't [fs]', \
                                ylab = 'E [V/m]', \
                        title = 'IFFT');

            ax3.plot(tt, IFT_drive, '-', label = 'Drive');
            ax3.plot(tt, IFT_wit, '-', label = 'Witness');
            ax3.plot(tt, IFT, '-.', label = 'Combined');
            ax3.legend();

            plt.show();
        return (FE_drive, FE_wit, FE), f;
    elif type(E) == np.ndarray or type(E) == list:
        FE, f, IFT, tt, f_plot = get_FFTs(E, t);
        
        if plot:
            NFFT = 2**(nextpow2(len(t)) + 1);
            fig1, ax1 = makefig(xlab = 't [fs]', \
                                ylab = 'E [V/m]', \
                        title = 'Input Field');
            ax1.plot(t * 1e15, E, label = 'Combined');
            
            fig2, ax2 = makefig(xlab = 'f [THz]', \
                                ylab = r'$F[E_r]$', \
                        title = 'FFT');
            FE_plot       = 2 * abs(FE[1:int(NFFT / 2 + 2)]);
            ax2.plot(f_plot, FE_plot, label = 'Combined');

            fig3, ax3 = makefig(xlab = 't [fs]', \
                                ylab = 'E [V/m]', \
                        title = 'IFFT');
            ax3.plot(tt, IFT, label = 'Combined');
            
            plt.show();
        return FE, f;
    else:
        print("ERROR: incorrect input electric field, returning null")
        return;
        

def cry_field(t, FEr, f, y0, ctype, d, nslice = 100, \
              plot = False, verbose = False):
    '''
    Calculates the THz field as it propagates through an EO crystal

    Parameters:
    -----------
    t      : array_like
             Time array for the electric field (s). 
    FEr    : array_like 
             Either an array of the FFT of the THz field or tuple of the FFTs  
             of the drive and witness fields
    f      : array_like
             Frequency array corresponding to FEr (THz);
    y0     : float
             The central wavelength of the probe pulse (m)
    ctype  : str
             The EO crystal being used (GaP or ZnTe)
    d      : float
             The crystal thickness (m)
    nslice : int, optional
             How many slices to divide the crystal into, default = 100
    plot   : bool, optional
             Whether or not to plot the field at certain points in the crystal
    Returns:
    --------
    '''

    if type(FEr) == tuple:
        N = 2;
        FE_drive = FEr[0]; 
        FE_wit   = FEr[1];
        FE       = FEr[2];
        # Only take forward propagating FFT
        sFEr_drive = FE_drive[0:int(len(FE) / 2 + 1)];
        sFEr_wit   = FE_wit[0:int(len(FE) / 2 + 1)];
        sFEr       = FE[0:int(len(FE) / 2 + 1)];
        sf   = f[0:int(len(FE) / 2 + 1)];
        f_Hz = sf * 1e12;
        # Compute transmission attenuation, pockel coeff, ind of ref
        eo_params     = {'ctype' : ctype, 'f' : f_Hz}
        A             = ec.transcoeff(eo_params);
        eps, n, kappa = ec.dielec(eo_params);
        r41           = ec.eocoeff(eo_params);
        
        
        # Propagate fields through crystal
        FEc_drive, Ec_drive = prop_field(sFEr_drive, \
                                         A, r41, f_Hz, n, kappa, d, t, nslice);
        FEc_wit, Ec_wit     = prop_field(sFEr_wit, \
                                         A, r41, f_Hz, n, kappa, d, t, nslice);
        FEc, Ec             = prop_field(sFEr, \
                                         A, r41, f_Hz, n, kappa, d, t, nslice);
        # Propagate reflection through crystal
        A_ref     = (1 - n - 1j * kappa) / (1 + n + 1j * kappa);
        
        F2_drive, E_ref_drive = prop_ref(sFEr_drive, A, r41,\
                                         A_ref, f_Hz, n, kappa, d, t, nslice);
        F2_wit, E_ref_wit     = prop_ref(sFEr_wit, A, r41, \
                                         A_ref, f_Hz, n, kappa, d, t, nslice);
        F2, E_ref             = prop_ref(sFEr, A, r41,\
                                         A_ref, f_Hz, n, kappa, d, t, nslice); 
        # Reconstruct time array
        Nsamp = len(sf);
        fmax  = sf[-1] * 1e12;
        tt    = np.arange(0, Nsamp - 1) * (1 / (fmax));
        
        tdiff = len(tt) - len(t);
        
        Ec_drive = Ec_drive[0:int(len(tt) - tdiff + 1), :];
        Ec_wit   = Ec_wit[0:int(len(tt) - tdiff + 1), :];
        Ec       = Ec[0:int(len(tt) - tdiff + 1), :];
        
        E_ref_drive = E_ref_drive[0:int(len(tt) - tdiff + 1), :];
        E_ref_wit   = E_ref_wit[0:int(len(tt) - tdiff + 1), :];
        E_ref       = E_ref[0:int(len(tt) - tdiff + 1), :];
        
        Ec_drive = Ec_drive# + E_ref_drive;
        Ec_wit   = Ec_wit#   + E_ref_wit;
        Ec       = Ec #+ E_ref;
        
        tt = tt[0:int(len(tt) - tdiff + 1)];
        tt = tt * 1e12; # ps
        #tt    = tt - 3 * tt[-1] / 4;
    elif type(FEr) == np.ndarray or type(FEr) == list:
        N = 1;
        FE     = FEr;
        # Only take forward propagating FFT
        sFEr       = FE[0:int(len(FE) / 2 + 1)];
        sf   = f[0:int(len(FEr) / 2 + 1)];
        f_Hz = sf * 1e12;
        # Compute transmission attenuation, pockel coeff, ind of ref
        eo_params     = {'ctype' : ctype, 'f' : f_Hz}
        A             = ec.transcoeff(eo_params);
        eps, n, kappa = ec.dielec(eo_params);
        r41           = ec.eocoeff(eo_params);
        

        
        
        # Propagate fields through crystal
        FEc, Ec             = prop_field(sFEr, \
                                         A, r41, f_Hz, n, kappa, d, t, nslice);
        
        # Propagate reflection through crystal
        A_ref     = (1 - n - 1j * kappa) / (1 + n + 1j * kappa);
        
        F2, E_ref             = prop_ref(sFEr, A, r41,\
                                         A_ref, f_Hz, n, kappa, d, t, nslice); 
        # Reconstruct time array
        Nsamp = len(sf);
        fmax  = sf[-1] * 1e12;
        tt    = np.arange(0, Nsamp - 1) * (1 / (fmax));
        tdiff = len(tt) - len(t);
        Ec       = Ec[0:int(len(tt) - tdiff + 1), :];
        E_ref    = E_ref[0:int(len(tt) - tdiff + 1), :];
        Ec       = Ec + E_ref;
        
        tt = tt[0:int(len(tt) - tdiff + 1)];
        tt = tt * 1e12; # ps
        tt    = tt - 3 * tt[-1] / 4;
    else:
        print("ERROR: incorrect input electric field, returning null")
        return; 
    
    # Create array of depths
    j      = np.arange(1, nslice, 1);
    dz     = d / nslice;
    d_arr  = (j - 0.5) * dz;
    # Print statements
    print_vel(ctype, f_Hz, verbose = verbose);
    t_shift, vg_eff = shift_probe(y0, Ec, tt, d, nslice, ctype, \
                              verbose = verbose);
    # Plot if requested:
    if plot:
        
        plot_probe = True;
        # Plotting all 3 gets very messy, just plot combined
        height = 20;
        fig1 =  plt.figure(figsize = (6, height), dpi = 200);
        
        # Get indices to plot (5 points in the crystal)
        z_step  = int(d * 1e6 / 5) # microns
        z_start = (z_step / 2);
        z_end   = (z_step * (5 - 0.5));
        z_plot  = np.arange(z_start, z_end + 1, z_step) * 1e-6;
        
        
        # Plot, bottom to top
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[4]));
        ax1           = fig1.add_subplot(515);
        ax1.set_prop_cycle(cy)
        Ec_plot       = np.real(Ec[:, ind]);
        Ec_wit_plot   = np.real(Ec_wit[:, ind])
        Ec_drive_plot = np.real(Ec_drive[:, ind])
        ax1.plot(tt, Ec_drive_plot, '-', label = 'drive')
        ax1.plot(tt, Ec_wit_plot, '-', label = 'wit')
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax1.plot(tt, Ec_plot, '-', label = 'combined');
        ax1.set_xlabel('t [ps]');
        # Peak location
        if plot_probe:
            t_probe = (z_plot[4] / vg_eff) * 1e12
            ax1.axvline(x = t_probe + t_shift, color = 'r');
            
        # Pretty plot stuff
        ax1.set_yticks([]);
        ax1.text(0.95, 0.85, depth, ha='center', va='center', \
                 transform=ax1.transAxes)
        
        ########################################################################
        ind = np.argmin(abs(d_arr - z_plot[3]));
        ax2 = fig1.add_subplot(514, sharex = ax1);
        ax2.set_prop_cycle(cy)
        Ec_plot       = np.real(Ec[:, ind]);
        Ec_wit_plot   = np.real(Ec_wit[:, ind])
        Ec_drive_plot = np.real(Ec_drive[:, ind])
        ax2.plot(tt, Ec_drive_plot, '-', label = 'drive')
        ax2.plot(tt, Ec_wit_plot, '-', label = 'wit')
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax2.plot(tt, Ec_plot, '-', label = r'$E_{eff}$');
        if plot_probe:
            t_probe = (z_plot[3] / vg_eff) * 1e12
            ax2.axvline(x = t_probe + t_shift, color = 'r');
        # Pretty plot stuff
        ax2.text(0.95, 0.85, depth, ha='center', va='center', \
                 transform=ax2.transAxes)
        ax2.set_yticks([]);
        
        
        ########################################################################
        ind = np.argmin(abs(d_arr - z_plot[2]));
        ax3 = fig1.add_subplot(513, sharex = ax1);
        ax3.set_prop_cycle(cy)
        Ec_plot       = np.real(Ec[:, ind]);
        Ec_wit_plot   = np.real(Ec_wit[:, ind])
        Ec_drive_plot = np.real(Ec_drive[:, ind])
        ax3.plot(tt, Ec_drive_plot, '-', label = 'drive')
        ax3.plot(tt, Ec_wit_plot, '-', label = 'wit')
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax3.plot(tt, Ec_plot, '-', label = r'$E_{eff}$');

        # Peak location
        if plot_probe:
            t_probe = (z_plot[2] / vg_eff) * 1e12
            ax3.axvline(x = t_probe + t_shift, color = 'r');
        # Pretty plot stuff
        ax3.text(0.95, 0.85, depth, ha='center', va='center', \
                 transform=ax3.transAxes)
        ax3.set_yticks([]);

        ########################################################################
        ind = np.argmin(abs(d_arr - z_plot[1]));
        ax4 = fig1.add_subplot(512, sharex = ax1);
        ax4.set_prop_cycle(cy)
        Ec_plot       = np.real(Ec[:, ind]);
        Ec_wit_plot   = np.real(Ec_wit[:, ind])
        Ec_drive_plot = np.real(Ec_drive[:, ind])
        ax4.plot(tt, Ec_drive_plot, '-', label = 'drive')
        ax4.plot(tt, Ec_wit_plot, '-', label = 'wit')
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax4.plot(tt, Ec_plot, '-', label = r'$E_{eff}$');
        # Peak location
        if plot_probe:
            t_probe = (z_plot[1] / vg_eff) * 1e12
            ax4.axvline(x = t_probe + t_shift, color = 'r');
        # Pretty plot stuff
        ax4.text(0.95, 0.85, depth, ha='center', va='center', \
                 transform=ax4.transAxes)
        ax4.set_yticks([]);

        ########################################################################
        ind = np.argmin(abs(d_arr - z_plot[0]));
        ax5 = fig1.add_subplot(511, sharex = ax1);
        ax5.set_prop_cycle(cy)
        Ec_plot       = np.real(Ec[:, ind]);
        Ec_wit_plot   = np.real(Ec_wit[:, ind])
        Ec_drive_plot = np.real(Ec_drive[:, ind])
        ax5.plot(tt, Ec_drive_plot, '-', label = 'drive')
        ax5.plot(tt, Ec_wit_plot, '-', label = 'wit')
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax5.plot(tt, Ec_plot, linestyle = '-',\
                 label = r'$E_{eff}$');
        
        ax5.legend(ncol = 1, loc = 'upper left', \
              bbox_to_anchor = (-0.15, 1.0), frameon = False);
        
        if plot_probe:
            t_probe = (z_plot[0] / vg_eff) * 1e12
            ax5.axvline(x = t_probe + t_shift, color = 'r', label = 'probe');
            
        # Pretty plot stuff
        ax5.set_yticks([]);
        #ind = np.argmax(Ec_wit_plot)
        ax5.text(0.95, 0.85, depth, ha='center', va='center', \
                 transform=ax5.transAxes);
        ########################################################################
        plt.setp(ax2.get_xticklabels(), visible=False);
        plt.setp(ax3.get_xticklabels(), visible=False);
        plt.setp(ax4.get_xticklabels(), visible=False);
        plt.setp(ax5.get_xticklabels(), visible=False);
        plt.subplots_adjust(hspace = 0.0);
        plt.show();
    if N == 1:
        return Ec, tt, FEc, sf, t_shift;
    elif N == 2:
        Ecs = (Ec_drive, Ec_wit, Ec);
        FEcs = (FEc_drive, FEc_wit, FEc);
        return Ecs, tt, FEcs, sf, t_shift;