import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
from scipy.special import gamma as gm
from scipy.integrate import simps
sys.path.insert(0, "../")
import Constants.SI as SI
global plasmaDict, c
c = SI.lightSpeed


class beam:
	def __init__(self, en = [5.3e-6], \
		sigma_z = [5.2e-6], sigma_r = [5.1e-6], Q = [1.5e-9], npoints = 100):
		''' Initializes an electron beam with Facet II like parameters
		Parameters:
		-----------
		gamma : array_like, optional
		Relativistic factor for the centroid beam energy, default 20000.
		en : array_like, optional
			Normalized emittance of the beam [mm-rad], default 5.3e-6
		beta_s : array_like, optional
			Array of beam waist betas [m], default [.1]
		sigma_z : array_like, optional
			Longitudinal beam size [m], default 5.2e-6
		Q : array_like, optional
			Charge of the beam [C], default 1.5e-9
	'''
		self.en = np.array(en);
		self.sigma_z = np.array(sigma_z);
		self.Q = np.array(Q);
		# Relativistic beta factor
		#self.beta = np.sqrt(1 - 1 / self.gamma**2)
		# beam sigma_r
		self.sigma_r = np.array(sigma_r)
		# (r,z) position of max E field
		self.r = np.zeros((npoints,1))
		self.z = np.zeros((1,npoints))
		self.r[:,0] = np.linspace(-50e-6, 50e-6, npoints)
		self.z[0,:] = np.linspace(-30e-6, 30e-6, npoints)
	def peak_charge_dens(self):
		self.ppk = self.Q / ((2*np.pi)**(3/2) * self.sigma_r**2 \
								* self.sigma_z) 

	def max_frac(self, Vi, Z,  eps0 = SI.permFreeSpace):
		self.max_ion_frac = np.zeros((len(self.sigma_r), len(self.sigma_z)))
		Vh = 13.6;
		n = Z / np.sqrt(Vi/Vh)
		t = self.z /  c;
		for i in range(len(self.sigma_r)):
			for j in range(len(self.sigma_z)):
				Er = np.squeeze(abs((self.ppk[i,j] * \
						self.sigma_r[i]**2 / (eps0 * self.r)) * \
						(1 - np.exp(-self.r**2/(2*self.sigma_r[i]**2))) *\
						np.exp(-(self.z)**2 / (2 * self.sigma_z[j]**2))))
				W  = 1.52 * ((4**n * Vi) / (n * gm(2*n))) * \
					 (20.5 * Vi**(3/2) / Er)**(2*n-1) * \
					 np.exp(-6.83 * Vi**1.5 /Er);
		
				plasma_frac = 1 - np.exp(-simps(W * 1e15, t));
				self.max_ion_frac[i,j] = np.max(plasma_frac);
		
