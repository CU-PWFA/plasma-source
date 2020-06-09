'''
Module for all things EO crystal
'''
import matplotlib.pyplot as plt;
import numpy as np;
from scipy.constants import c;
from scipy.constants import epsilon_0 as eps0;

from plotting import makefig;
class crystal:

    def __init__(self, ctype):
        # ctype : str, crystal type either GaP or ZnTe

        self.ctype = ctype;

        if self.ctype.lower() == 'gap':
            self.epsel  = 8.7;
            self.f0     = 10.98e12; # Hz
            self.s0     = 1.8;
            self.gamma0 = 0.02e12; # Hz
            self.dE     = 1e-12; 
            self.C      = -0.53;
        elif self.ctype.lower() == 'znte':
            self.epsel  = 7.4;
            self.f0     = 5.3e12;
            self.s0     = 2.7;
            self.gamma0 = 0.09e12; 
            self.dE     = 4.25e-12;
            self.C      = -0.07;
        else:
            print('Crystal type not specified')
            self.epsel  = 1;
            self.f0     = 0;
            self.s0     = 0;
            self.gamma0 = 0;
            self.dE     = 1;
            self.C      = 0;
    def dielec(self, f, plot = False, log = False):
        ''' Calculates the complex dielectric function for a range of THz
            frequencies 
    
        Parameters:
        -----------
        f     : array_like
                The range of frequencies to compute the dielectric 
                function
                in Hz
        plot  : bool, optional
                Whether or not to plot results
        log   : bool, optional
                Whether or not to have the xscale of plots be log
        '''

        eps   = self.epsel + self.s0 * self.f0**2 \
                / (self.f0**2 - f**2 - 1j * self.gamma0 * f); # complex dielec
        n     = np.real(np.sqrt(eps)); # index of refraction
        kappa = np.imag(np.sqrt(eps)); # attenuation coefficient;


        if plot:
            fig1, ax1 = makefig(xlab = 'f [THz]', ylab = 'n(f)', \
                                title = self.ctype + ' index of refraction');
            ax1.plot(f * 1e-12, n)

            fig2, ax2 = makefig(xlab = 'f [THz]', ylab = r'$\kappa$(f)', \
                                title = self.ctype + ' attenuation');
            ax2.plot(f * 1e-12, kappa)

            plt.show();
        return eps, n, kappa

    def transcoeff(self, f, plot = False):
        '''
        Function to calculate the transmission coefficient as a function of 
        frequency. 
    
        Parameters:
        -----------
        f     : array_like
                The range of frequencies to compute the dielectric function
                in Hz
        plot  : bool, optional
                Whether or not to plot results
   
        Returns:
        --------
        A : array_like
            The transmission coefficient 
        '''
        eps, n, kappa = self.dielec(f);
        A = 2 / (1 + n + 1j * kappa);

        if plot:
            fig, ax = makefig(xlab = 'f [THz]', ylab = r'$A_{tr}$', \
                              title = self.ctype + ' transmission coefficient');
            ax.plot(f * 1e12, A)
            plt.show();
        return A
    def eocoeff(self, f, plot = False):
        '''
        Function to compute the pockels coefficient of the crystal as a function
        of frequency. 
        Parameters:
        -----------
        f : array_like
            Frequency in Hz
        plot : bool, optional
               Whether or not to plot the pockel coefficient
        Returns:
        --------
        r41 : array_like
              The pockel coefficient
        '''

        r41 = self.dE * (1 + self.C * self.f0**2 \
              / (self.f0**2 - f**2 - 1j * self.gamma0 * f));

        if plot:
            fig, ax = makefig(xlab = 'f [THz]', ylab = r'$r_{41}$ [m/v]', \
                      title = self.ctype + ' pockel coefficient');
            ax.plot(f * 1e-12, r41)
            plt.show()
        return r41
    def indref(self, y_arr, plot = False):
        '''
        Function to calculate the optical index of refrection of a crystal as a 
        function of wavelength. 
        
        Parameters:
        -----------
        ctype : str
                Specifies the EO crystal, either GaP or ZnTe
        y_arr : array_like
                Wavelength in m
        plot  : bool, optional
                Whether or not to plot results
        Returns:
        --------
        n : array_like
            index of refraction
        '''


        # convert y_arr to microns
        lam = y_arr * 1e6;

        if self.ctype.lower() == 'gap':
            n = np.sqrt(2.68 + 6.40 * lam**2 / (lam**2 - 0.0903279));
        elif self.ctype.lower() == 'znte':
            und = np.argwhere(lam < 30);
            ove = np.argwhere(lam >= 30);
            n = np.zeros(len(lam));
            n[und] = np.sqrt(9.92 + 0.42530 / (lam[und]**2 - 0.37766**2) \
                         + 2.63580 / ((lam[und]**2) / (56.5**2) - 1));
            n[ove] = np.sqrt(4.270 + 3.01 * lam[ove]**2 / (lam[ove]**2 -0.142));

        else:
            n = np.zeros(len(lam)) + 1.0;
        if plot:
            fig, ax = makefig(xlab = r'$\lambda$ [nm]', ylab = 'n', \
                              title = self.ctype + ' index of refraction');
            ax.plot(lam * 1e3, n);
            plt.show()
        return n;

    def velocities(self, f, plot = False):
        '''
        Function to calculate the phase velocity, group velocity, and group 
        velocity dispersion in the crystal as a function of frequency.

        Parameters:
        -----------
        f    : array_like
               Frequency in Hz
        plot : bool, optional
               Whether or not to plot results

        Returns:
        --------
        v_ph : array_like
               Phase velocity
        v_g  : array_like
               Group velocity
        g_vd : array_like
               Group velocity dispersion
        '''

        just_one = False;
        if len(f) == 1:
            just_one = True;
            f = np.array([0.999 * f[0], f[0], 1.001 * f[0]]);
        # Get index of refraction
        if all(i < 1e14 for i in f): # THz range
            eps, n, kappa = self.dielec(f);
        else:
            y_arr = (c / f); # m
            y_arr[f == 0] = np.nan;
            n = self.indref(y_arr);
        # get dn/df
        dndf_calc = np.diff(n) / np.diff(f);
        dndf      = np.append(dndf_calc, dndf_calc[-1]);

        if just_one:
            n = n[1];
            f = f[1];
            dndf = dndf[1];
        v_ph = c / n;
        v_g  = c / (n + f * dndf);
        g_vd = (2 / (c * np.pi)) * dndf * 1e15; # fs^2 / mm

        if plot:
            fig, ax = makefig(xlab = 'f [THz]', ylab = r'$\frac{v}{c}$', 
                              title = self.ctype + ' velocities');

            ax.plot(f * 1e-12, v_ph / c, label = r'$v_{\phi}$');
            ax.plot(f * 1e-12, v_g / c, label = r'$v_g$');

            # plot 800 nm group velocity for comparison
            y_opt   = np.array([0.999 * 800, 800, 1.001 * 800]) * 1e-9; # m
            f_opt     = c / y_opt;
            n_opt     = self.indref(y_opt);
            dndf_calc = np.diff(n_opt) / np.diff(f_opt);
            dndf      = np.append(dndf_calc, dndf_calc[-1]);
            n_opt     = n_opt[1];
            f_opt     = f_opt[1];
            v_g_opt   = c / (n_opt + f_opt * dndf[1]);
            ax.axhline(v_g_opt / c, linestyle = '-.', color = 'k', \
                       label = r'$v_g$ (800 nm)');
            ax.legend();
            ax.set_ylim([0, v_ph[0] / c])
            plt.show();
        return v_ph, v_g, g_vd

    def georesponse(self, f, d, probe, plot = False):
        '''
        Function to compute the geometric response function of the crystal as a 
        function of frequency and crystal thickness

        Parameters:
        -----------
        f     : array_like
                Frequency in Hz
        d     : array_like
                Crystal thicknesses in m
        probe : object
                An instance of the probe class
        plot  : bool, optional
                Whether or not to plot the geometric response function

        Returns:
        G : array_like
            The crystal geometric response function
        '''

        # THz range
        v_ph, v_g, g_vd = self.velocities(f);
        
        # Optical
        lam_opt = np.linspace(probe.y0 - 25e-9, probe.y0 + 25e-9, 51); # m
        f_opt   = (c / lam_opt); # Hz
        v_ph_opt, v_g_opt, d_vg_opt = self.velocities(f_opt);
        v_g_opt = np.mean(v_g_opt);
        
        # Attenuation coeff
        eps, n, kappa = self.dielec(f);
        a = 2 * np.pi * f * kappa / c;
        # Set longitudinal distance array
        nz = 1000;
        G  = np.zeros((len(f), len(d)), dtype = 'complex128');
        
        inv_del_v = (1 / v_ph - 1/v_g_opt);
        
        for i in range(len(d)):
            z  = np.linspace(0, d[i], nz);
            dz = abs(z[1]-z[0]);
            for j in range(len(z)):
                G[:, i] += (1 / d[i]) \
                           * np.exp(1j * 2 * np.pi * f * z[j] \
                                    * inv_del_v) * np.exp(-a * z[j]) * dz;
        if plot:
            fig, ax = makefig(xlab = 'f [THz]', ylab = '|G(f)|', \
                              title = self.ctype + ' geometric response');
            for i in range(len(d)):
                lab = 'd = ' + str(np.round(d[i] * 1e6)) + r' $\mu$m';
                ax.plot(f * 1e-12, np.abs(G[:, i]), label = lab);
            ax.legend();
            plt.show();
        return G

    def eoresponse(self, f, d, probe, plot = False):
        '''
        Function to compute the geometric response function of the crystal as a 
        function of frequency and crystal thickness

        Parameters:
        -----------
        f     : array_like
                Frequency in Hz
        d     : array_like
                Crystal thicknesses in m
        probe : object
                An instance of the probe class
        plot  : bool, optional
                Whether or not to plot the geometric response function

        Returns:
        G_eo : array_like
               The crystal electro-optic response function
        '''

        G = self.georesponse(f, d, probe);

        A = self.transcoeff(f);

        r41 = self.eocoeff(f);

        G_eo = np.zeros((len(f), len(d)));

        if plot:
            fig, ax = makefig(xlab = 'f [THz]', ylab = r'$G_{EO}$ [pV / m]', \
                              title = self.ctype + ' electro-optic response')
        for i in range(len(d)):
            G_eo[:, i] = G[:, i] * A * r41;
            if plot:
                lab = 'd = ' + str(np.round(d[i] * 1e6)) + r' $\mu$m';
                ax.plot(f * 1e-12, np.abs(G_eo[:, i]) * 1e12, label = lab);
        if plot:
            ax.legend();
            plt.show();
        return G_eo;
