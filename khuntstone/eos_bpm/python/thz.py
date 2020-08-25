'''
Module for propagating a THz pulse through an EO crystal
'''
from cycler import cycler;
import matplotlib.pyplot as plt;
import numpy as np;
from scipy.constants import c;
from scipy.constants import epsilon_0 as eps0;
import scipy.fftpack as sfft;

from plotting import makefig;

# colors for plotting
plot_colors = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', \
               '#DDCC77', '#CC6677', '#882255', '#AA4499'];
cy = cycler('color', plot_colors);
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

    L  = len(t);
    Ts = abs(t[1] - t[0]);
    Fs = 1 / Ts;
    def nextpow2(i):
        n = 1;
        while 2**n < i : n += 1
        return n

    # Get the FFT
    NFFT   = 2**(nextpow2(L) + 1);
    f      = Fs * np.linspace(0, 1, NFFT) * 1e-12; # THz frequencies
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
    tt    = tt * 1e12; # ps
    
    if plot:
            NFFT = 2**(nextpow2(len(t)) + 1);
            fig1, ax1 = makefig(xlab = 't [ps]', \
                                ylab = 'E [V/m]', \
                        title = 'Input Field');
            ax1.plot(t * 1e12, E, label = 'Combined');
            
            fig2, ax2 = makefig(xlab = 'f [THz]', \
                                ylab = r'$F[E_r]$', \
                        title = 'FFT');
            FE_plot       = 2 * abs(FEr[1:int(NFFT / 2 + 2)]);
            ax2.plot(f_plot, FE_plot, label = 'Combined');

            fig3, ax3 = makefig(xlab = 't [ps]', \
                                ylab = 'E [V/m]', \
                        title = 'IFFT');
            ax3.plot(tt, IFER, label = 'Combined');
            
            plt.show();
    return FEr, f;

def cry_field(t, FEr, f, d, probe, crystal, nslice = 100, plot = False, \
              verbose = False, save = False, sname = ''):
    '''
    Function to calculate the THz electric field at each slice in the crystal

    Parameters:
    -----------
    t       : array_like
              Time array for the electric field (s)
    FEr     : array_like
              FFT of the THz pulse
    f       : array_like
              Frequency array corresponding to FEr (Hz)
    d       : float
              Crystal thickness (m)
    probe   : object
              Instance of laser class
    nslice  : int, optional
              Number of slices to divide the crystal into
    plot    : bool, optional
              Whether or not to plot the propagation 
    verbose : bool, optional
              Whether or not to surpress print statements.
    save    : bool, optional
              Whether or not to save the plot
    sname   : str, optional
              The name to save the figure
    '''

    # Only take forward propagating FFT

    sFEr = FEr[0:int(len(FEr) / 2 + 1)];
    sf   = f[0:int(len(FEr) / 2 + 1)];

    f_Hz = sf * 1e12; 

    # crystal parameters
    A             = crystal.transcoeff(f_Hz);
    eps, n, kappa = crystal.dielec(f_Hz);
    r41           = crystal.eocoeff(f_Hz);
    #r41           = 1;

    # Divide crystal into slices
    j     = np.arange(1, nslice, 1);
    dz    = d / nslice;
    d_arr = (j - 0.5) * dz;

    # Preallocate for loop
    FEc = np.zeros((len(f_Hz), len(d_arr)), dtype = 'complex128');

    Ec  = np.zeros((len(f_Hz), len(d_arr)), dtype = 'complex128');

    # Wavenumber
    k = 2 * np.pi * f_Hz * n / c;
    # Attenuation
    alpha = 2 * np.pi * f_Hz * kappa / c

    # Propagate pulse

    for i in range(len(d_arr)):
        FEc[:, i] = A * sFEr * r41 * np.exp((1j * 2 * np.pi * f_Hz / c) \
                      * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                      * d_arr[i]);

        Ec[:, i] = sfft.ifft(FEc[:, i]) * len(t);
        Ec[:, i] = np.flip(Ec[:, i], axis = 0);
    

    A_ref     = (1 - n - 1j * kappa) / (1 + n + 1j * kappa);
    k         = 2 * np.pi * f_Hz * n / c; 
    alpha     = 2 * np.pi * f_Hz * kappa / c;
    ref_phase = np.exp(1j * 2 * k * d) * np.exp(-2 * alpha * d);
    F2       = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
    Eref       = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');

    for i in range(len(d_arr)):           
        F2[:, i]       = A * sFEr * r41 * A_ref**2 * ref_phase  \
                           * np.exp((1j * 2 * np.pi * f_Hz / c) \
                           * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                           * d_arr[i]);                          
        Eref[:, i]       = sfft.ifft(F2[:, i]) * len(t);
        Eref[:, i]       = np.flip(Eref[:, i], axis = 0);

    Nsamp = len(sf); 
    fmax  = sf[-1] * 1e12;
    tt    = np.arange(0, Nsamp - 1) * (1 / (fmax));
    tdiff = len(tt) - len(t);

    Ec    = Ec[0:int(len(tt) - tdiff + 1), :];
    Eref  = Eref[0:int(len(tt) - tdiff + 1), :];
    
    Ec    = Ec + Eref;
    tt    = tt[0:int(len(tt) - tdiff + 1)];
    tt    = tt * 1e12;
    tt    = tt - tt[int( 3 * len(tt) / 4) + 1];
    
    
    # Add time shift so probe is matched to THz peak in the first slice
    y_opt   = np.array([0.999, 1, 1.001]) * probe.y0;
    f_opt     = c / y_opt;
    n_opt     = crystal.indref(y_opt);
    dndf_calc = np.diff(n_opt) / np.diff(f_opt);
    dndf      = np.append(dndf_calc, dndf_calc[-1]);
    n_opt     = n_opt[1];
    f_opt     = f_opt[1];
    v_g_opt   = c / (n_opt + f_opt * dndf[1]);
    match_ind = 0;
    t_peak = tt[np.argmax(np.real(Ec[:, match_ind]))];
    probe.t_shift = t_peak - ((d_arr[match_ind] / v_g_opt) * 1e12);
    probe.t_shift = 0
    #print(probe.t_shift)
    if verbose:
        v_ph, v_g, dummy = crystal.velocities(f_Hz);
        print("THz phase velocity:", np.round(np.nanmean(v_ph) / c, 2), 'c');
        print("THz group velocity:", np.round(np.nanmean(v_g) / c, 2), 'c');
        print("Probe group velocity:", np.round(np.nanmean(v_g_opt) / c, 2), 'c');
        print("Probe shift:", np.round(probe.t_shift, 2), 'ps')

    if plot:
        # Get characterstic broadening length
        Lchar, dummy = probe.laser_l_char(crystal);
        height = 20;
        fig    = plt.figure(figsize = (8, height), dpi = 200);

        # get indices to plot (5 points in crystal)

        z_step = int(d * 1e6 / 5); # microns
        z_start = z_step / 2;
        z_end   = z_step * (4.5);
        z_plot  = np.arange(z_start, z_end + 1, z_step) * 1e-6;

        # Plot, bottom to top
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[4]));
        ax1           = fig.add_subplot(515);
        ax1.set_prop_cycle(cy)
        Ec_plot       = np.real(Ec[:, ind]) / max(np.real(Ec[:, ind]));
        depth         = f'{np.round(d_arr[ind] * 1e6, 0)}' + r'$\mu$m';
        ax1.plot(tt, Ec_plot, '-', label = depth);
        # Plot the probe
        t_plot_probe = (d_arr[ind] / v_g_opt) * 1e12 + probe.t_shift;
        sigj         = probe.sigp * np.sqrt(1 + (d_arr[ind] / Lchar)**2);

        t_exp        = (tt - t_plot_probe) * 1e-12
        probe_plot   = np.exp(-(t_exp)**2 / (2 * sigj**2));
        ax1.plot(tt, probe_plot, '-r')
       
        ax1.set_xlabel('t [ps]');
        

        lg1 = ax1.legend(handlelength=0, handletextpad=0, fancybox=True);
        for item in lg1.legendHandles:
            item.set_visible(False)
        # Pretty plot stuff
        ax1.set_yticks([]);
        ind2 = np.argmax(np.real(Eref[:, ind]));
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[3]));
        ax2           = fig.add_subplot(514);
        ax2.set_prop_cycle(cy)
        Ec_plot       = np.real(Ec[:, ind]) / max(np.real(Ec[:, ind]));
        depth         = f'{np.round(d_arr[ind] * 1e6, 0)}' + r'$\mu$m';
        ax2.plot(tt, Ec_plot, '-', label = depth);

        # Plot the probe
        t_plot_probe = (d_arr[ind] / v_g_opt) * 1e12 + probe.t_shift;
        t_exp        = (tt - t_plot_probe) * 1e-12
        probe_plot   = np.exp(-(t_exp)**2 / (2 * sigj**2));
        ax2.plot(tt, probe_plot, '-r')


        lg2 = ax2.legend(handlelength=0, handletextpad=0, fancybox=True);
        for item in lg2.legendHandles:
            item.set_visible(False)
        # Pretty plot stuff
        ax2.set_yticks([]);
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[2]));
        ax3           = fig.add_subplot(513);
        ax3.set_prop_cycle(cy)
        Ec_plot       = np.real(Ec[:, ind]) / max(np.real(Ec[:, ind]));
        depth         = f'{np.round(d_arr[ind] * 1e6, 0)}' + r'$\mu$m';
        ax3.plot(tt, Ec_plot, '-', label = depth);

        # Plot the probe
        t_plot_probe = (d_arr[ind] / v_g_opt) * 1e12 + probe.t_shift;
        t_exp        = (tt - t_plot_probe) * 1e-12
        probe_plot   = np.exp(-(t_exp)**2 / (2 * sigj**2));
        ax3.plot(tt, probe_plot, '-r')


        lg3 = ax3.legend(handlelength=0, handletextpad=0, fancybox=True);
        for item in lg3.legendHandles:
            item.set_visible(False)
        # Pretty plot stuff
        ax3.set_yticks([]);
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[1]));
        ax4           = fig.add_subplot(512);
        ax4.set_prop_cycle(cy)
        Ec_plot       = np.real(Ec[:, ind]) / max(np.real(Ec[:, ind]));
        depth         = f'{np.round(d_arr[ind] * 1e6, 0)}' + r'$\mu$m';
        ax4.plot(tt, Ec_plot, '-', label = depth);

        # Plot the probe
        t_plot_probe = (d_arr[ind] / v_g_opt) * 1e12 + probe.t_shift;
        t_exp        = (tt - t_plot_probe) * 1e-12
        probe_plot   = np.exp(-(t_exp)**2 / (2 * sigj**2));
        ax4.plot(tt, probe_plot, '-r')


        lg4 = ax4.legend(handlelength=0, handletextpad=0, fancybox=True);
        for item in lg4.legendHandles:
            item.set_visible(False)
        # Pretty plot stuff
        ax4.set_yticks([]);
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[0]));
        ax5           = fig.add_subplot(511);
        ax5.set_prop_cycle(cy)
        Ec_plot       = np.real(Ec[:, ind]) / max(np.real(Ec[:, ind]));
        depth         = f'{np.round(d_arr[ind] * 1e6, 0)}' + r'$\mu$m';
        ax5.plot(tt, Ec_plot, '-', label = depth);

        # Plot the probe
        t_plot_probe = (d_arr[ind] / v_g_opt) * 1e12 + probe.t_shift;
        t_exp        = (tt - t_plot_probe) * 1e-12
        probe_plot   = np.exp(-(t_exp)**2 / (2 * sigj**2));
        ax5.plot(tt, probe_plot, '-r')

        
        lg5 = ax5.legend(handlelength=0, handletextpad=0, fancybox=True);
        for item in lg5.legendHandles:
            item.set_visible(False)
        ax5.set_yticks([]);
        ind1 = np.argmax(Ec_plot);
        ########################################################################


        # Prettify
        t_min = tt[0]
        t_max = tt[-1]
        
        ax1.set_xlim([t_min, t_max]);
        ax2.set_xlim([t_min, t_max]);
        ax3.set_xlim([t_min, t_max]);
        ax4.set_xlim([t_min, t_max]);
        ax5.set_xlim([t_min, t_max]);

        ax5.set_xticks([])
        ax4.set_xticks([])
        ax3.set_xticks([])
        ax2.set_xticks([])

        plt.subplots_adjust(hspace = 0)
        fig.savefig(sname);
    return Ec, tt;


def cry_two_field(t, FEr_drive, FEr_wit, FEr, f, d, probe, crystal,\
                  nslice = 100, plot = False, verbose = False):
    '''
    Similar function to cry_field, but now computes the individual fields of a 
    drive and witness beam

    Parameters:
    -----------
    t         : array_like
                Time array for the electric field (s)
    FEr_drive : array_like
                FFT of the drive THz pulse
    FEr_wit   : array_like
                FFT of the witness THz pulse
    FEr       : array_like
                FFT of the combined fields
    f         : array_like
                Frequency array corresponding to FEr (Hz)
    d         : float
                Crystal thickness (m)
    probe     : object
                Instance of laser class
    nslice    : int, optional
                Number of slices to divide the crystal into
    plot      : bool, optional
                Whether or not to plot the propagation 
    verbose   : bool, optional
                Whether or not to surpress print statements.
    '''

    # Only take forward propagating FFT

    sFEr_drive = FEr_drive[0:int(len(FEr) / 2 + 1)];
    sFEr_wit   = FEr_wit[0:int(len(FEr) / 2 + 1)];
    sFEr       = FEr[0:int(len(FEr) / 2 + 1)];
    sf         = f[0:int(len(FEr) / 2 + 1)];

    f_Hz = sf * 1e12; 

    # crystal parameters
    A             = crystal.transcoeff(f_Hz);
    eps, n, kappa = crystal.dielec(f_Hz);
    r41           = crystal.eocoeff(f_Hz);


    # Divide crystal into slices
    j     = np.arange(1, nslice, 1);
    dz    = d / nslice;
    d_arr = (j - 0.5) * dz;

    # Preallocate for loop
    FEc_drive = np.zeros((len(f_Hz), len(d_arr)), dtype = 'complex128');
    FEc_wit   = np.zeros((len(f_Hz), len(d_arr)), dtype = 'complex128');
    FEc       = np.zeros((len(f_Hz), len(d_arr)), dtype = 'complex128');


    Ec_drive  = np.zeros((len(f_Hz), len(d_arr)), dtype = 'complex128');
    Ec_wit    = np.zeros((len(f_Hz), len(d_arr)), dtype = 'complex128');
    Ec        = np.zeros((len(f_Hz), len(d_arr)), dtype = 'complex128');

    # Wavenumber
    k = 2 * np.pi * f_Hz * n / c;
    # Attenuation
    alpha = 2 * np.pi * f_Hz * kappa / c

    # Propagate pulse

    for i in range(len(d_arr)):
        FEc_drive[:, i] = A * sFEr_drive * r41 \
                          *  np.exp((1j * 2 * np.pi * f_Hz / c) \
                          * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                          * d_arr[i]);

        FEc_wit[:, i]   = A * sFEr_wit * r41 \
                          * np.exp((1j * 2 * np.pi * f_Hz / c) \
                          * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                          * d_arr[i]);

        FEc[:, i]       = A * sFEr * r41 * np.exp((1j * 2 * np.pi * f_Hz / c) \
                          * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                          * d_arr[i]);

        Ec_drive[:, i] = sfft.ifft(FEc_drive[:, i]) * len(t);
        Ec_drive[:, i] = np.flip(Ec_drive[:, i], axis = 0);

        Ec_wit[:, i] = sfft.ifft(FEc_wit[:, i]) * len(t);
        Ec_wit[:, i] = np.flip(Ec_wit[:, i], axis = 0);

        Ec[:, i] = sfft.ifft(FEc[:, i]) * len(t);
        Ec[:, i] = np.flip(Ec[:, i], axis = 0);
    

    A_ref     = (1 - n - 1j * kappa) / (1 + n + 1j * kappa);
    k         = 2 * np.pi * f_Hz * n / c; 
    alpha     = 2 * np.pi * f_Hz * kappa / c;
    ref_phase = np.exp(1j * 2 * k * d) * np.exp(-2 * alpha * d);


    F2_drive   = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
    Eref_drive = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');

    F2_wit     = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
    Eref_wit   = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');

    F2         = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');
    Eref       = np.zeros((len(sf), len(d_arr)), dtype = 'complex128');

    for i in range(len(d_arr)):  
        F2_drive[:, i] = A * sFEr_drive * r41 * A_ref**2 * ref_phase  \
                           * np.exp((1j * 2 * np.pi * f_Hz / c) \
                           * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                           * d_arr[i]);  

        F2_wit[:, i]   = A * sFEr_wit * r41 * A_ref**2 * ref_phase  \
                           * np.exp((1j * 2 * np.pi * f_Hz / c) \
                           * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                           * d_arr[i]);  

        F2[:, i]       = A * sFEr * r41 * A_ref**2 * ref_phase  \
                           * np.exp((1j * 2 * np.pi * f_Hz / c) \
                           * n * d_arr[i] - (2 * np.pi * f_Hz / c) * kappa \
                           * d_arr[i]); 
        Eref_drive[:, i]       = sfft.ifft(F2_drive[:, i]) * len(t);
        Eref_drive[:, i]       = np.flip(Eref_drive[:, i], axis = 0);

        Eref_wit[:, i]       = sfft.ifft(F2_wit[:, i]) * len(t);
        Eref_wit[:, i]       = np.flip(Eref_wit[:, i], axis = 0);  

        Eref[:, i]       = sfft.ifft(F2[:, i]) * len(t);
        Eref[:, i]       = np.flip(Eref[:, i], axis = 0);

    Nsamp = len(sf); 
    fmax  = sf[-1] * 1e12;
    tt    = np.arange(0, Nsamp - 1) * (1 / (fmax));
    tdiff = len(tt) - len(t);

    Ec_drive    = Ec_drive[0:int(len(tt) - tdiff + 1), :];
    Eref_drive  = Eref_drive[0:int(len(tt) - tdiff + 1), :];
    Ec_drive    = Ec_drive + Eref_drive;

    Ec_wit      = Ec_wit[0:int(len(tt) - tdiff + 1), :];
    Eref_wit    = Eref_wit[0:int(len(tt) - tdiff + 1), :];
    Ec_wit      = Ec_wit + Eref_wit;

    Ec          = Ec[0:int(len(tt) - tdiff + 1), :];
    Eref        = Eref[0:int(len(tt) - tdiff + 1), :];
    Ec    = Ec + Eref;



    # Reconstruct time array
    tt    = tt[0:int(len(tt) - tdiff + 1)];
    tt    = tt * 1e12;
    tt    = tt - tt[int( 3 * len(tt) / 4) + 1];
    
    
    # Add time shift so probe is matched to THz peak in the first slice
    y_opt   = np.array([0.999, 1, 1.001]) * probe.y0;
    f_opt     = c / y_opt;
    n_opt     = crystal.indref(y_opt);
    dndf_calc = np.diff(n_opt) / np.diff(f_opt);
    dndf      = np.append(dndf_calc, dndf_calc[-1]);
    n_opt     = n_opt[1];
    f_opt     = f_opt[1];
    v_g_opt   = c / (n_opt + f_opt * dndf[1]);
    match_ind = 0;
    t_peak = tt[np.argmax(np.real(Ec[:, match_ind]))];
    probe.t_shift = t_peak - ((d_arr[match_ind] / v_g_opt) * 1e12);

    if verbose:
        v_ph, v_g, dummy = crystal.velocities(f_Hz);
        print("THz phase velocity:", np.round(np.mean(v_ph) / c, 2), 'c');
        print("THz group velocity:", np.round(np.mean(v_g) / c, 2), 'c');
        print("Probe group velocity:", np.round(np.mean(v_g_opt) / c, 2), 'c');
        print("Probe shift:", np.round(probe.t_shift, 2), 'ps')

    if plot:
        # Get characterstic broadening length
        Lchar, dummy = probe.laser_l_char(crystal);
        height = 20;
        fig    = plt.figure(figsize = (8, height), dpi = 200);

        # get indices to plot (5 points in crystal)

        z_step = int(d * 1e6 / 5); # microns
        z_start = z_step / 2;
        z_end   = z_step * (4.5);
        z_plot  = np.arange(z_start, z_end + 1, z_step) * 1e-6;

        # Plot, bottom to top
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[4]));
        ax1           = fig.add_subplot(515);
        ax1.set_prop_cycle(cy)
        Ec_drive_plot = np.real(Ec_drive[:, ind]) / max(np.real(Ec[:, ind]));
        Ec_wit_plot   = np.real(Ec_wit[:, ind]) / max(np.real(Ec[:, ind]));
        Ec_plot       = np.real(Ec[:, ind]) / max(np.real(Ec[:, ind]));
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax1.plot(50, 50, '-w', label = 'd = ' + depth);
        ax1.plot(tt, Ec_plot, '-', label = 'Combined');
        ax1.plot(tt, Ec_drive_plot, '--', label = 'Drive');
        ax1.plot(tt, Ec_wit_plot, '-.', label = 'Witness');

        # Plot the probe
        t_plot_probe = (d_arr[ind] / v_g_opt) * 1e12 + probe.t_shift;
        sigj         = probe.sigp * np.sqrt(1 + (d_arr[ind] / Lchar)**2);

        t_exp        = (tt - t_plot_probe) * 1e-12
        probe_plot   = np.exp(-(t_exp)**2 / (2 * sigj**2));
        ax1.plot(tt, probe_plot * max(Ec_plot), '-r')
        
        ax1.set_xlabel('t [ps]');
        

        ax1.legend();
        # Pretty plot stuff
        ax1.set_yticks([]);
        ind2 = np.argmax(Ec_plot);
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[3]));
        ax2           = fig.add_subplot(514);
        ax2.set_prop_cycle(cy)
        Ec_drive_plot = np.real(Ec_drive[:, ind]) / max(np.real(Ec[:, ind]));
        Ec_wit_plot   = np.real(Ec_wit[:, ind]) / max(np.real(Ec[:, ind]));
        Ec_plot       = np.real(Ec[:, ind]) / max(np.real(Ec[:, ind]));
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax2.plot(50, 50, '-w', label = 'd = ' + depth);
        ax2.plot(tt, Ec_plot, '-', label = 'Combined');
        ax2.plot(tt, Ec_drive_plot, '--', label = 'Drive');
        ax2.plot(tt, Ec_wit_plot, '-.', label = 'Witness');

        # Plot the probe
        t_plot_probe = (d_arr[ind] / v_g_opt) * 1e12 + probe.t_shift;
        t_exp        = (tt - t_plot_probe) * 1e-12
        probe_plot   = np.exp(-(t_exp)**2 / (2 * sigj**2));
        ax2.plot(tt, probe_plot * max(Ec_plot), '-r')

        ax2.set_xlabel('t [ps]');
        ax2.legend();
        # Pretty plot stuff
        ax2.set_yticks([]);
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[2]));
        ax3           = fig.add_subplot(513);
        ax3.set_prop_cycle(cy);
        Ec_drive_plot = np.real(Ec_drive[:, ind]) / max(np.real(Ec[:, ind]));
        Ec_wit_plot   = np.real(Ec_wit[:, ind]) / max(np.real(Ec[:, ind]));
        Ec_plot       = np.real(Ec[:, ind]) / max(np.real(Ec[:, ind]));
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax3.plot(50, 50, '-w', label = 'd = ' + depth);
        ax3.plot(tt, Ec_plot, '-', label = 'Combined');
        ax3.plot(tt, Ec_drive_plot, '--', label = 'Drive');
        ax3.plot(tt, Ec_wit_plot, '-.', label = 'Witness');

        # Plot the probe
        t_plot_probe = (d_arr[ind] / v_g_opt) * 1e12 + probe.t_shift;
        t_exp        = (tt - t_plot_probe) * 1e-12
        probe_plot   = np.exp(-(t_exp)**2 / (2 * sigj**2));
        ax3.plot(tt, probe_plot * max(Ec_plot), '-r')

        ax3.set_xlabel('t [ps]');
        ax3.legend();
        # Pretty plot stuff
        ax3.set_yticks([]);
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[1]));
        ax4           = fig.add_subplot(512);
        ax4.set_prop_cycle(cy)
        Ec_drive_plot = np.real(Ec_drive[:, ind]) / max(np.real(Ec[:, ind]));
        Ec_wit_plot   = np.real(Ec_wit[:, ind]) / max(np.real(Ec[:, ind]));
        Ec_plot       = np.real(Ec[:, ind]) / max(np.real(Ec[:, ind]));
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax4.plot(50, 50, '-w', label = 'd = ' + depth);
        ax4.plot(tt, Ec_plot, '-', label = 'Combined');
        ax4.plot(tt, Ec_drive_plot, '--', label = 'Drive');
        ax4.plot(tt, Ec_wit_plot, '-.', label = 'Witness');

        # Plot the probe
        t_plot_probe = (d_arr[ind] / v_g_opt) * 1e12 + probe.t_shift;
        t_exp        = (tt - t_plot_probe) * 1e-12
        probe_plot   = np.exp(-(t_exp)**2 / (2 * sigj**2));
        ax4.plot(tt, probe_plot * max(Ec_plot), '-r')

        ax4.set_xlabel('t [ps]');
        ax4.legend();
        # Pretty plot stuff
        ax4.set_yticks([]);
        ########################################################################
        ind           = np.argmin(abs(d_arr - z_plot[0]));
        ax5           = fig.add_subplot(511);
        ax5.set_prop_cycle(cy)
        Ec_drive_plot = np.real(Ec_drive[:, ind]) / max(np.real(Ec[:, ind]));
        Ec_wit_plot   = np.real(Ec_wit[:, ind]) / max(np.real(Ec[:, ind]));
        Ec_plot       = np.real(Ec[:, ind]) / max(np.real(Ec[:, ind]));
        depth         = f'{np.round(d_arr[ind] * 1e6, 2)}' + r'$\mu$m';
        ax5.plot(50, 50, '-w', label = 'd = ' + depth);
        ax5.plot(tt, Ec_plot, '-', label = 'Combined');
        ax5.plot(tt, Ec_drive_plot, '--', label = 'Drive');
        ax5.plot(tt, Ec_wit_plot, '-.', label = 'Witness');

        # Plot the probe
        t_plot_probe = (d_arr[ind] / v_g_opt) * 1e12 + probe.t_shift;
        t_exp        = (tt - t_plot_probe) * 1e-12
        probe_plot   = np.exp(-(t_exp)**2 / (2 * sigj**2));
        ax5.plot(tt, probe_plot * max(Ec_plot), '-r')

        ax5.set_xlabel('t [ps]');
        ax5.legend();
        ax5.set_yticks([]);
        Eref_plot = np.real(Eref[:, ind])
        ind1 = np.argmax(Ec_plot);
        ########################################################################


        # Prettify
        t_min = tt[ind1] - 0.5;
        t_max = tt[ind2] + 0.5;


        y_min = 1.2 * np.amin(np.real(Ec_plot));
        y_max = 1.2 * np.amax(np.real(Ec_plot));
        ax1.set_ylim([y_min, y_max]);
        ax2.set_ylim([y_min, y_max]);
        ax3.set_ylim([y_min, y_max]);
        ax4.set_ylim([y_min, y_max]);
        ax5.set_ylim([y_min, y_max]);

        ax1.set_xlim([t_min, t_max]);
        ax2.set_xlim([t_min, t_max]);
        ax3.set_xlim([t_min, t_max]);
        ax4.set_xlim([t_min, t_max]);
        ax5.set_xlim([t_min, t_max]);

        ax5.set_xticks([])
        ax4.set_xticks([])
        ax3.set_xticks([])
        ax2.set_xticks([])
    return Ec_drive, Ec_wit, Ec, tt;
