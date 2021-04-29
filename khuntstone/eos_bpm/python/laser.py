'''
Class for all things probe laser.

author: keenan
'''

import matplotlib.pyplot as plt;
import numpy as np;
from scipy.constants import c;
from scipy.interpolate import interp1d;
from plotting import makefig;


class laser:
    # Object for all probe pulse parameters and functions

    keys = ['y0', 'tp', 'dy'];

    def __init__(self, params):
        '''
        Initializes probe class

        Parameters:
        -----------
        params : dictionary
                 dictionary containing keys 'y0', 'tp', and 'dy' with values
                 for central wavelength (m), FWHM (s), and bandwidth (m) 
                 respectively
        '''

        self.params = params;
        self.check_params(params);
        self.params_to_attr(params);

        # Get central frequency and frequency bandwidth
        self.w0 = 2 * np.pi * c / self.y0;
        y1      = self.y0 - self.dy / 2;
        y2      = self.y0 + self.dy / 2;
        w1      = 2 * np.pi * c / y1;
        w2      = 2 * np.pi * c / y2;
        self.dw = abs(w1 - w2);
        # Get rms
        self.sigp = self.tp / (2 * np.sqrt(2 * np.log(2)));

    def check_params(self, params):
        for key in self.keys:
            if key not in params:
                raise ValueError('The params object has no key %s' % key);
    def params_to_attr(self, params):
        '''
        Add all params as attributes of the class
        '''
        for key in params:
            setattr(self, key, params[key]);

    def chirp(self, tc):
        '''
        Chirps the probe pulse
        Parameters:
        -----------
        tc : float
             The FWHM duration of the chirped pulse (s)
        '''

        self.chirped = True;
        self.tc      = tc;

        # Compute chirp parameters a, b and Gamma:

        self.a  = 2 * np.log(2) / (self.tc**2);
        self.b  = np.sqrt(((self.a * self.dw**2) \
                     / (8 * np.log(2))) - self.a**2);
        self.gm = self.a + 1j * self.b;

        
    def check_chirp(self):
        '''
        Checks if pulse has been chirped.
        '''
        if hasattr(self, 'chirped'):
            return True;
        else:
            return False
    def get_inst_w(self, t):
        '''
        Computes the instantaneous frequency of the chirped pulse at time t
        Parameters:
        -----------
        t : float or array_like
            The time at which to compute the instantaneous frequency (s) 
            

        Returns:
        --------
        wi : float or array_like 
             The instantaneous frequency
        '''
        if not self.check_chirp():
            raise ValueError('Pulse is not chirped');
            return;
        else:
            return self.w0  + 2 * self.b * t;
    def get_inst_amp(self, t):
        '''
        Computes the instantaneous amplitude (AU, max 1) of the chirped pulse
        Parameters:
        -----------
        t : float or array_like
            The time at which to compute the instantaneous frequency (s)

        Returns:
        --------
        A : float or array_like 
             The instantaneous amplitude
        '''
        if not self.check_chirp():
            raise ValueError('Pulse is not chirped');
            return;
        else:
            t_arr  = np.linspace(-10 * self.tc, 10 * self.tc, 1000);
            Ec_env = np.exp(-self.gm * t_arr**2) *\
                            np.exp(-1j * self.w0 * t_arr);
            Ec_env = np.abs(Ec_env)**2;
            f      = interp1d(t_arr, Ec_env);
            A      = f(t)
            return A

    def plot_field(self):
        if self.check_chirp():
            fig, ax = makefig(xlab = 't [ps]', ylab = r'|E|$^2$ [AU]', title = 
                  'Probe field envelope');
            t_plot = np.linspace(-self.tc, self.tc, 1000);
            sigp   = self.tp / (2 * np.sqrt(2 * np.log(2)));
            E_env  = np.exp(-t_plot**2 / (2 * sigp **2));

            Ec_env = np.exp(-self.gm * t_plot**2) *\
                     np.exp(-1j * self.w0 * t_plot);
            Ec_env = np.abs(Ec_env)**2;
            t_plot = t_plot * 1e12;
            ax.plot(t_plot, E_env, label = r'$E_p$');
            ax.plot(t_plot, Ec_env, label = r'$E_c$');
            ax.legend();
            plt.show();
        else:
            fig, ax = makefig(xlab = 't [fs]', ylab = r'|E|$^2$ [AU]', title = 
                  'Probe field envelope');
            t_plot = np.linspace(-self.tp, self.tp, 1000);
            sigp   = self.tp / (2 * np.sqrt(2 * np.log(2)));
            E_env  = np.exp(-t_plot**2 / (2 * sigp **2));
            t_plot = t_plot * 1e15;
            ax.plot(t_plot, E_env);
            plt.show();
    def plot_minis(self):
        '''
        Function to plot the chirped pulse envelope and the mini-Gaussians it 
        is divided into. 
        '''
        if self.check_chirp():
            t_plot = np.linspace(-self.tc, self.tc, 1000);
            Ec_env = np.exp(-self.gm * t_plot**2) \
            * np.exp(-1j * self.w0 * t_plot);
            Ai    = self.get_inst_amp(t_plot);
            wi    = self.get_inst_w(t_plot); 
            mini_sig = np.sqrt(self.tp * self.tc)
            fig, ax = makefig(xlab = 't [ps]', ylab = r'$|E|^2$')
            for i in range(len(t_plot)):
                if i%10==0:
                    mini = Ai[i] \
                    * np.exp(-(t_plot - t_plot[i])**2 / (2 * mini_sig**2));
                    ax.plot(t_plot * 1e12, mini, '--r')
            ax.plot(t_plot * 1e12, np.abs(Ec_env)**2, linewidth = 8)
            plt.show()
        else:
            print("Probe pulse is unchirped");
    def laser_l_char(self, crystal):
        '''
        Function to compute the characteristic pulse broadening length of the 
        pulse in an EO crystal.

        Parameters:
        -----------
        Crystal : object
                  Instance of the crystal class
        Returns:
        --------
        L_char : float
                 The characteristic pulse broadening length (m)
        '''
        y_arr     = np.array([0.99, 1, 1.01]) * self.y0;
        f_opt     = c / y_arr;
        n_opt     = crystal.indref(y_arr)
        dndf_calc = np.diff(n_opt) / np.diff(f_opt);
        dndf      = np.append(dndf_calc, dndf_calc[-1]);
        vg_opt    = c / (n_opt + f_opt * dndf)
        w_opt     = 2 * np.pi * f_opt;
        dvg_dw    = np.diff(vg_opt) / np.diff(w_opt);
        L_char    = -2 * self.sigp**2 * vg_opt[1]**2 * (dvg_dw[0])**(-1);
        return L_char, vg_opt[1];


