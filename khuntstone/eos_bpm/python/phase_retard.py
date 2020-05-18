'''
Module to compute the phase retardation experienced by a probing pulse 
propagating through an EO crystal. Allows for several different probe setups

author : keenan
'''
import matplotlib.pyplot as plt;
import numpy as np;
from scipy.constants import c;
from scipy.interpolate import interp1d;

from laser import laser;
from plotting import makefig;
def phase_retard(Ec, te, d, tau, probe, crystal, method,\
                 psi = 0, interp = False, n_interp = 10000, plot = False):
    '''
    Function to compute phase retardation experienced by a probe pulse 
    propagating through an electo-optic crystal

    Parameters:
    -----------
    Ec : array_like
         The temporal profile of the THz electric field in each slice of the 
         crstal
    te : array_like
         The time array corresponding to Ec (s)
    d  : float
            The crystal thickness (m)
    tau   : array_like
            Either the probe delay time (if ultrashort pulse) or the time array
            corresponding to the chirped probe pulse
    probe : object
            instance of the probe class
    crystal : object
              instance of the crystal class (be sure to use the same crystal
              that was used in calculating Ec)
    method : str
             The method for treating the probe pulse either 'broad' to 
             incorporate pulse broadening, 'chirp' for a chirped pulse,
             'delta' for a delta like pulse, and 'spatial' for spatial encoding
    psi     : float, optional
              The angle between the probe and crystal in degs, only used if 
              method is spatial
    interp  : bool, optional
              Whether or not to interpolate gamma
    n_interp : int, optional
               number of points for interpolating gamma (tau should be course
               to cut down on run time, final phase shift is interpolated)

    plot   : bool, optional
             Whether or not to plot the phase retardation

    Returns:
    --------
    gamma : array_like
            The phase retardation
    tg    : array_like
            The time array corresponding to gamma
    '''

    # Convert psi to rads
    psi = psi * np.pi / 180;
    # Create depth array;
    nslice = np.shape(Ec)[1] + 1;
    j     = np.arange(1, nslice, 1);
    dz    = d / nslice;
    d_arr = (j - 0.5) * dz;
    if method == 'broad':
        n0 = crystal.indref(np.array([probe.y0]))[0]
        dz = abs(d_arr[1] - d_arr[0])
        amp = 2 * np.pi * n0**3 * dz / (probe.y0 * np.sqrt(2 * np.pi));
        Lchar, vg_opt = probe.laser_l_char(crystal);
        gamma = np.zeros(len(tau));
        t_probe = probe.t_shift * 1e-12; # s
        for i in range(len(tau)):
            # Updates are nice
            if (i+1) % 100 == 0:
                print(np.round(i / len(tau) * 100), "%")
            for j in range(len(d_arr)):
                sigj = probe.sigp * np.sqrt(1 + (d_arr[j] / Lchar)**2)
                f    = interp1d(te, np.real(Ec[:, j]),\
                                bounds_error = False, fill_value = 0);
               
                t_interp = te + (d_arr[j] / vg_opt);
                E_eff = f(t_interp);
                t_delay = tau[i] + t_probe
                integrand = amp * E_eff * \
                            np.exp(-(te - t_delay)**2 / (2 * sigj**2)) / (sigj);
                gamma[i] += np.trapz(integrand, te);
        tau_interp = np.linspace(tau[0], tau[-1], 1000);   
        f          = interp1d(tau, gamma);     
        gamma_interp = f(tau_interp)
    
    elif method == 'chirp':
        t_probe = probe.t_shift * 1e-12;
        gamma  = np.zeros(len(tau));
        for i in range(len(tau)):
            # Updates are nice
            if (i+1) % 100 == 0:
                print(np.round(i / len(tau) * 100), "%")
            # Get mini-gaussian
            wi = probe.get_inst_w(tau[i]);
            yi = 2 * np.pi * c / wi
            sig0 = np.sqrt(probe.tc * probe.tp)
            mini_probe = laser({'y0' : yi, 'dy': 0,
                                'tp' : 2 * np.sqrt(2 * np.log(2)) * sig0});
            n0 = crystal.indref(np.array([mini_probe.y0]))[0];
            dz = abs(d_arr[0] - d_arr[1])
            amp = 2 * np.pi * n0**3 * dz / (mini_probe.y0 * np.sqrt(2 * np.pi));
            Lchar, vg_opt = mini_probe.laser_l_char(crystal);
            A = probe.get_inst_amp(tau[i]);
            for j in range(len(d_arr)):
                sigj = mini_probe.sigp * np.sqrt(1 + (d_arr[j] / Lchar)**2)
                f    = interp1d(te, np.real(Ec[:, j]),\
                                bounds_error = False, fill_value = 0);
               
                t_interp = te + (d_arr[j] / vg_opt);
                E_eff = f(t_interp);
                t_delay = tau[i] + t_probe
                integrand = amp * E_eff * \
                            np.exp(-(te - t_delay)**2 / (2 * sigj**2)) / (sigj);
                gamma[i] += np.trapz(integrand, te);
        tau_interp = np.linspace(tau[0], tau[-1], n_interp);
        f          = interp1d(tau, gamma);     
        gamma_interp = f(tau_interp)

    elif method == 'delta':
        n0 = crystal.indref(np.array([probe.y0]))[0]
        dz = abs(d_arr[1] - d_arr[0])
        amp = 2 * np.pi * n0**3 * dz / (probe.y0);
        Lchar, vg_opt = probe.laser_l_char(crystal);
        gamma = np.zeros(len(tau));
        t_probe = probe.t_shift * 1e-12; # s
        if psi != 0:
            t_delay = tau;
            vg_eff = vg_opt * np.cos(psi);
        else:
            t_delay = tau #+ t_probe;
            vg_eff = vg_opt;
        vg_eff  = vg_opt * np.cos(psi);
        for j in range(len(d_arr)):
            f    = interp1d(te, np.real(Ec[:, j]),\
                            bounds_error = False, fill_value = 0);
            t_interp = (d_arr[j] / vg_eff) + t_delay;
            E_eff = f(t_interp);
            gamma += amp * E_eff;

        if psi == 0:
            tau_interp = np.linspace(tau[0], tau[-1], 1000);   
            f          = interp1d(tau, gamma);     
            gamma_interp = f(tau_interp)
        else:
            tau_interp = tau;
            gamma_interp = gamma;

    elif method == 'spatial':
        n0 = crystal.indref(np.array([probe.y0]))[0];
        dz = abs(d_arr[1] - d_arr[0])
        amp = 2 * np.pi * n0**3 * dz / (probe.y0);
        gamma = np.zeros(len(tau));
        Lchar, v_g_opt = probe.laser_l_char(crystal);
        #print(max(tau), probe.t_shift*1e-12)
        tau = tau + probe.t_shift*1e-12
        vg_eff  = v_g_opt * np.cos(psi);
        for j in range(len(d_arr)):
            f = interp1d(te, np.real(Ec[:, j]), \
                        bounds_error = False, fill_value = 0);
            t_interp = tau + (d_arr[j] / vg_eff);
            E_eff = f(t_interp);
            gamma += amp * E_eff;
        tau_interp = np.linspace(tau[0], tau[-1], n_interp);
        f_gamma = interp1d(tau, gamma);
        gamma_interp = f_gamma(tau_interp);

    if plot:
        fig, ax = makefig(xlab = 't [fs]', ylab = r'$\Gamma$', \
                          title = 'Phase retardation')
        ax.plot(tau_interp * 1e15, gamma_interp);
        plt.show();
    if interp:
        return gamma_interp, tau_interp;
    else:
        return gamma, tau