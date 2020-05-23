'''
Module for developing current_profiles.py, testing functions, etc... 
'''

# Standard imports
from cycler import cycler
import matplotlib.pyplot as plt;
import numpy as np;
from scipy.signal import find_peaks, savgol_filter;
from scipy.interpolate import interp1d as interp
import sys;
sys.path.insert(0, "../../python/");
# Custom modules
import current_profiles as cp; # For reading in Claudio's data
import eosim as sim; # For running EOS-BPM simulations
import thz_prop as prop;
plot_colors = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', \
               '#DDCC77', '#CC6677', '#882255', '#AA4499'];

cy = cycler('color', plot_colors);

def makefig(x = 4, y = 3, xlab = '', ylab = '', title = '', fs = 12, ts = 12):
    fig = plt.figure(figsize = (x, y), dpi = 200);
    ax  = fig.gca();
    ax.set_prop_cycle(cy)
    ax.set_xlabel(xlab, fontsize = fs);
    ax.set_ylabel(ylab, fontsize = fs);
    ax.set_title(title, fontsize = ts, fontweight = 'bold');
    ax.tick_params(labelsize = 'large');
    return fig, ax;

ind = 0;
r   = 5e-3;
E, ze, te, I, ti = cp.get_E(ind, r);    
E_smooth = cp.smooth_E(E);  
sigt = 17e-15;
E_drive = 1e9 * np.exp(-te**2 / (2 * sigt**2)); 
#E_drive = np.zeros(len(te));
E_wit   = 0.3e9 * np.exp(-(te - 584e-15)**2 / (2 * sigt**2));         
#E_drive, E_wit = cp.split_E(E_test, te); 
#E_wit = np.zeros(len(te))                    
E_sim = (E_drive, E_wit, E_drive + E_wit);
tau   = np.linspace(-200, 1800, 1000) * 1e-15;
probe = {'y0' : 800e-9, 'width' : 27e-9, 'a_laser' : 0};
ctype = 'ZnTe'
d     = 100e-6;
ts    = 0;
gamma, t_plot = sim.sim(E_sim, te, tau, ctype, d, probe, shift = ts, \
                    x = [r, r], y = [0, 0], plot = True, plot_input = False, \
                    normed = True);
