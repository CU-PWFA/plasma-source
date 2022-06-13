'''
Electron beam class for computing THz electric field (assumed Gaussian) from
beam parameters.

author : keenan
'''
import numpy as np;
from scipy.constants import c;
from scipy.constants import epsilon_0 as eps0;
class ebeam:
    keys = ['Q', 'sigz', 'del_z', 't']

    def __init__(self, params):
        '''
        Initializes ebeam class

        Parameters:
        -----------
        params : dictionary
                 dictionary containing keys 'Q', 'sigz', 'del_z', t, and 'r' 
                 with values for charge (C), rms length (m),
                 longitudinal offset (0 for drive beam, m), time domain of the 
                 beam (s).
        '''

        self.params = params;
        self.check_params(params);
        self.params_to_attr(params);


        # Get temporal rms length and longitudinal offset
        self.sigt  = self.sigz / c;
        self.del_t = self.del_z / c;
        
    def get_Er(self, r):
        
        '''
        Computes the radial electric field of the bunch
        Params:
        -------
        r : float
            The radial distance from the electron bunch
        '''
        self.r = r;
        # Compute THz electric field of the beam
        self.E0 = self.Q / (2 * np.pi * eps0 * self.r * np.sqrt(2 * np.pi) \
         * c * self.sigt);

        self.Er = self.E0 * np.exp(-(self.t + self.del_t)**2 \
                                  / (2 * self.sigt**2));



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

