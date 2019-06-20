#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 10:49:43 2017

@author: robert
"""

import numpy as np
from beam.beams import beam
from vsim import load as Vload
from vsim import analyze as Vanalyze
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy import optimize
import scipy.integrate as Int

class ElectronBeam(beam.Beam):
    """ An electron beam class that stores a collection of macroparticles.
    
    Parameters
    ----------
    N : int
        The number of macroparticles in the beam.
    """
    keys = ['N']
    
    # Initialization functions
    #--------------------------------------------------------------------------
    
    def __init__(self, params):
        super().__init__(params)
        # Create internal variables
        if self.load is False:
            self.initialize_particles()
            self.save_initial()
    
    def initialize_particles(self, ptcls=None):
        """ Create the arrays to store particle positions and momentum in. 
        
        Parameters
        ----------
        ptcls : array-like, optional
            The array of particle position momentum to initialize.
        """
        if ptcls is None:
            self.ptcls = np.zeros((self.N, 6), dtype='double')
        else:
            self.ptcls = np.array(ptcls, dtype='double')
        self.saveInd = 0
        self.z = []
        self.save_ptcls(self.ptcls, 0.0)
    
    def load_beam(self):
        """ Load the beam, specfically load the z-grid and saveInd. """
        self.z = np.load(self.filePre + '_z.npy')
        self.saveInd = len(self.z)
        self.ptcls = self.load_ptcls(self.saveInd - 1)[0]
    
    # Getters and setters
    #--------------------------------------------------------------------------
    
    def get_x(self, ptcls):
        return ptcls[:, 0]
    
    def get_xp(self, ptcls):
        return ptcls[:, 1]
    
    def get_y(self, ptcls):
        return ptcls[:, 2]
    
    def get_yp(self, ptcls):
        return ptcls[:, 3]
    
    def get_z(self, ptcls):
        return ptcls[:, 4]
    
    def get_gamma(self, ptcls):
        return ptcls[:, 5]
    
    #File managment
    #--------------------------------------------------------------------------
    
    def save_ptcls(self, ptcls, z):
        """ Save the current particle distribution to file and advance z. """
        np.save(self.filePre + '_ptcls_' + str(self.saveInd) + '.npy', ptcls)
        self.saveInd += 1
        self.z.append(z)
        self.save_z()
    
    def save_z(self):
        """ Save the z array. """
        np.save(self.filePre + '_z.npy', self.z)
    
    def load_ptcls(self, ind):
        """ Load the particle distribution at the specified index. 
        
        Parameters
        ----------
        ind : int
            The save index to load the field at.
        
        Returns
        -------
        ptcls : array-like
            The particle corrdinates and momenta at the specified index.
        z : double
            The z coordinate of the field.
        """
        ptcls = np.load(self.filePre + '_ptcls_' + str(ind) + '.npy')
        z = self.z[ind]
        return ptcls, z
    
    def get_save_z(self,ind):
        return self.load_ptcls(ind)[1]
    
    # Physics functions
    #--------------------------------------------------------------------------
        
    def propagate(self, z, n):
        """ Propagate the field to an array of z distances. """
        #TODO implement this function
    
    def get_emittance(self, ind, ptcls=None, weights=None):
        """ Calculate the emittance from a particular save file. """
        ptcls = self.load_ptcls(ind)[0]
        x = self.get_x(ptcls)
        xp = self.get_xp(ptcls)
        y = self.get_y(ptcls)
        yp = self.get_yp(ptcls)
        #If weights aren't given, just initialize an array of 1's
        if weights is None:
            weights = np.zeros(len(x))+1
        # Calculate the differences from the average
        dx = x - np.average(x, weights=weights)
        dxp = xp - np.average(xp, weights=weights)
        dy = y - np.average(y, weights=weights)
        dyp = yp - np.average(yp, weights=weights)
        # Calculate the RMS sizes and the correlation
        sigmax2 = np.average(dx**2, weights=weights)
        sigmaxp2 = np.average(dxp**2, weights=weights)
        sigmaxxp = np.average(dx*dxp, weights=weights)
        sigmay2 = np.average(dy**2, weights=weights)
        sigmayp2 = np.average(dyp**2, weights=weights)
        sigmayyp = np.average(dy*dyp, weights=weights)
        # Calculate the emittance
        ex = np.sqrt(sigmax2*sigmaxp2 - sigmaxxp**2)
        ey = np.sqrt(sigmay2*sigmayp2 - sigmayyp**2)
        return ex, ey
    
    def get_gamma_n(self, ind, weights=None):
        """ Calculate the Lorentz factor from a particular save file. """
        ptcls = self.load_ptcls(ind)[0]
        gamma_arr = self.get_gamma(ptcls)
        #If weights aren't given, just initialize an array of 1's
        if weights is None:
            weights = np.zeros(len(gamma_arr))+1
        gamma = np.average(gamma_arr, weights=weights)
        return gamma
    
    def get_emittance_n(self, ind, weights=None):
        """ Calculate the normalized emittance from a particular save file. """
        ex, ey = self.get_emittance(ind, weights=weights)
        gamma = self.get_gamma_n(ind, weights=weights)
        ex = ex*gamma
        ey = ey*gamma
        return ex, ey
    
    def get_sigmar(self, ind, weights=None):
        """ Caclulate the transverse beam size from a particular save file. """
        ptcls = self.load_ptcls(ind)[0]
        x = self.get_x(ptcls)
        y = self.get_y(ptcls)
        #If weights aren't given, just initialize an array of 1's
        if weights is None:
            weights = np.zeros(len(x))+1
        # Calculate the differences from the average
        dx = x - np.average(x, weights=weights)
        dy = y - np.average(y, weights=weights)
        # Calculate the RMS sizes and the correlation
        sigmax = np.sqrt(np.average(dx**2, weights=weights))
        sigmay = np.sqrt(np.average(dy**2, weights=weights))
        return sigmax, sigmay
    
    def get_sigmarp(self, ind, weights=None):
        """ Calculate the beam divergence from a particular save file. """
        ptcls = self.load_ptcls(ind)[0]
        xp = self.get_xp(ptcls)
        yp = self.get_yp(ptcls)
        #If weights aren't given, just initialize an array of 1's
        if weights is None:
            weights = np.zeros(len(xp))+1
        # Calculate the differences from the average
        dxp = xp - np.average(xp, weights=weights)
        dyp = yp - np.average(yp, weights=weights)
        # Calculate the RMS sizes and the correlation
        sigmaxp = np.sqrt(np.average(dxp**2, weights=weights))
        sigmayp = np.sqrt(np.average(dyp**2, weights=weights))
        return sigmaxp, sigmayp

    def get_beam_properties(self, ind):
        """ Calculate most of the beam properties from a save file. 
        
        Creates an output dictionary exactly matching Mike's code.
        """
        prop = {}
        
        ptcls = self.load_ptcls(ind)[0]
        x = self.get_x(ptcls)
        xp = self.get_xp(ptcls)
        gamma = np.average(self.get_gamma(ptcls))
        # Calculate the differences from the average
        dx = x - np.average(x)
        dxp = xp - np.average(xp)
        # Calculate the RMS sizes and the correlation
        sigmax2 = np.average(dx**2)
        sigmaxp2 = np.average(dxp**2)
        sigmaxxp = np.average(dx*dxp)
        # Calculate the emittance
        exn = gamma*np.sqrt(sigmax2*sigmaxp2 - sigmaxxp**2)
        prop['x_eps']   = exn
        prop['x']       = np.sqrt(sigmax2)
        prop['xp']      = np.sqrt(sigmaxp2)
        prop['xxp']     = np.sqrt(sigmaxxp)
        prop['x_beta']  = gamma*sigmax2/exn
        prop['x_gamma'] = gamma*sigmaxp2/exn
        prop['x_alpha'] = -gamma*sigmaxxp/exn
        prop['x_phase'] = np.arctan2(2*prop['x_alpha'], 
                          prop['x_gamma']-prop['x_beta'])/2
        return prop
    
    def get_CS_at(self, ind, weights=None):
        """ Find the CS parameters and beam parameters at a given index.
        
        Returns
        -------
        beam : dictionary
            eps_x, eps_y : double
                Geometric emittance in x and y.
            beta_x, beta_y : double
                Beta function in x and y.
            alpha_x, alpha_y : double
                Alpha function in x and y.
            gamma_x, gamma_y : double
                Gamma function in x and y.
            gamma_b : double
                Relativistic gamma of the beam.
            cen_x, cen_y : double
                The transverse beam center in x and y. 
        """
        beam = {}
        
        ptcls = self.load_ptcls(ind)[0]
        x = self.get_x(ptcls)
        xp = self.get_xp(ptcls)
        y = self.get_y(ptcls)
        yp = self.get_yp(ptcls)
        gamma = np.average(self.get_gamma(ptcls), weights=weights)
        # Calculate the differences from the average
        cen_x = np.average(x, weights=weights)
        dx = x - cen_x
        dxp = xp - np.average(xp, weights=weights)
        cen_y = np.average(y, weights=weights)
        dy = y - cen_y
        dyp = yp - np.average(yp, weights=weights)
        # Calculate the RMS sizes and the correlation
        sigmax2 = np.average(dx**2, weights=weights)
        sigmaxp2 = np.average(dxp**2, weights=weights)
        sigmaxxp = np.average(dx*dxp, weights=weights)
        sigmay2 = np.average(dy**2, weights=weights)
        sigmayp2 = np.average(dyp**2, weights=weights)
        sigmayyp = np.average(dy*dyp, weights=weights)
        # Calculate the emittance
        ex = np.sqrt(sigmax2*sigmaxp2 - sigmaxxp**2)
        ey = np.sqrt(sigmay2*sigmayp2 - sigmayyp**2)
        beam['eps_x']   = ex
        beam['eps_y']   = ey
        beam['beta_x']  = sigmax2/ex
        beam['beta_y']  = sigmay2/ey
        beam['alpha_x'] = sigmaxxp/ex
        beam['alpha_y'] = sigmayyp/ey
        beam['gamma_x'] = sigmaxp2/ex
        beam['gamma_y'] = sigmayp2/ey
        beam['gamma_b'] = gamma
        beam['cen_x']   = cen_x
        beam['cen_y']   = cen_y
        return beam
    
    def get_CS(self, weights=None):
        """ Return arrays of the CS parameters in each direction.
        
        Returns
        -------
        beam : dictionary
            eps_x, eps_y : array of double
                Geometric emittance in x and y.
            beta_x, beta_y : array of double
                Beta function in x and y.
            alpha_x, alpha_y : array of double
                Alpha function in x and y.
            gamma_x, gamma_y : array of double
                Gamma function in x and y.
            gamma_b : array of double
                Relativistic gamma of the beam.
            cen_x, cen_y : double
                The transverse beam center in x and y. 
        """
        z = self.z
        N = len(z)
        beam = {}
        keys = ['eps_x', 'eps_y', 'beta_x', 'beta_y', 'alpha_x', 'alpha_y',
                'gamma_x', 'gamma_y', 'gamma_b', 'cen_x', 'cen_y']
        for key in keys:
            beam[key] = np.zeros(N, dtype='double')
        for i in range(N):
            step = self.get_CS_at(i, weights)
            for key, item in step.items():
                beam[key][i] = item
        return beam
    
    def get_ptcl(self, ind):
        """ Load the 6D phase space for a single particle in the beam.
        
        Parameters
        ----------
        ind : int
            The index of the beam particle.
        
        Returns
        -------
        ptcl : array of double
            An array with 6 rows describing the particles full phase space.
        """
        z = self.z
        N = len(z)
        ptcl = np.zeros((6, N), dtype='double')
        for i in range(N):
            step = self.load_ptcls(i)[0][ind, :]
            ptcl[:, i] = step[:6]
        return ptcl
        
    # Visualization functions
    #--------------------------------------------------------------------------
    
    def plot_current_phase(self, xlim=None, ylim=None):
        """ Plots a scatterplot of the particles in the current beam. """
        self.plot_phase(self.ptcls, self.z[-1], xlim, ylim)
        plt.show()
    
    def plot_phase_at(self, ind):
        """ Plots the particles at a particular z distance.
        
        Parameters
        ----------
        ind : int
            The save index to plot the particles at, see the _z file to find z.
        """
        ptcls, z = self.load_ptcls(ind)
        self.plot_phase(ptcls, z)
        plt.show()
    
    def plot_phase(self, ptcls, z, xlim=None, ylim=None, weights = None):
        """ Create an x-y plot of the particles. """        
        #If weights aren't given, just initialize an array of 1's
        if weights is None:
            weights = np.zeros(self.N)+1
        else:
            sort = np.argsort(weights)
            weights = weights[sort]
            ptcls = ptcls[sort]
            
        #zf = self.get_z(ptcls)
        #weights = zf
            
        fig = plt.figure(figsize=(10, 4), dpi=150)
        plt.subplot(121)
        sctx = plt.scatter(ptcls[:, 0]*1e6, ptcls[:, 1]*1e3, 1, c=weights)
        plt.title('x phase space')
        plt.xlabel('x (um)')
        plt.ylabel("x' (mrad)")
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)
        plt.subplot(122)
        scty = plt.scatter(ptcls[:, 2]*1e6, ptcls[:, 3]*1e3, 1, c=weights)
        plt.title('y phase space')
        plt.xlabel('y (um)')
        plt.ylabel("y' (mrad)")
        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)

        plt.tight_layout()
        return fig, sctx, scty
    
    def plot_hist_at(self, ind):
        ptcls, z = self.load_ptcls(ind)
        self.plot_hist(ptcls, z, ind)
    
    def plot_hist(self, ptcls, z, ind, weights = None):
        if weights is None:
            weights = np.zeros(self.N)+1
        else:
            sort = np.argsort(weights)
            weights = weights[sort]
            ptcls = ptcls[sort]
        
        sigma = self.get_sigmar(ind)[0]
        #There used to be many different fits, but now only the Gauss+Gauss remains.
        # If you want to implement more, use the general procudure below with
        # functions at the bottom of this script.
        
        plt.figure(figsize=(10, 4), dpi=150)
        plt.subplot(121)
        left = plt.hist(ptcls[:,0]*1e6, bins = 50, weights = weights, log=True)        
        x = np.linspace(min(ptcls[:,0]),max(ptcls[:,0]),100)
        fx = max(left[0])*np.exp(-1*np.square(x)/2/np.square(sigma))
        plt.plot(x*1e6, fx, label=r'$\sigma_x=$'+str(sigma*1e6)+r'$\ \mu m$')
        plt.xlabel(r'$x\mathrm{\ [\mu m]}$')
        plt.ylabel('Counts')
        plt.ylim(bottom = 0.1)
        plt.ylim(top = max(left[0])*6)
        plt.legend()
        
        plt.subplot(122)
        right = plt.hist(ptcls[:,0]*1e6, bins = 50, weights = weights, log=False)
        plt.plot(x*1e6, fx, label=r'$\sigma_x=$'+str(sigma*1e6)+r'$\ \mu m$')
        plt.xlabel(r'$x\mathrm{\ [\mu m]}$')
        plt.ylim(bottom = 0.1)
        plt.ylim(top = max(right[0])*1.2)
        plt.legend(); plt.show()
        
        #With the data from left we now want to try a Gauss + Gauss fit
        wdata = np.array(left[0])
        rdata = np.array(left[1])*1e-6
        rdata = rdata[0:-1]+.5*(rdata[1]-rdata[0])
        p = FitDataSomething(wdata, rdata, GaussPlusGauss, [sigma, sigma*3, max(left[0]), max(left[0])/100])
        print("G+G: ",p)
        
        plt.figure(figsize=(10, 4), dpi=150)
        plt.subplot(121)
        left = plt.hist(ptcls[:,0]*1e6, bins = 50, weights = weights, log=True)
        fxg = GaussPlusGauss(p,x)
        plt.plot(x*1e6, fxg, label="G+G: "+r'$\sigma_1=$'+("%0.2f"%(np.abs(p[0])*1e6))+r'$\ \mu m$'+" & "+r'$\sigma_2=$'+("%0.2f"%(np.abs(p[1])*1e6))+r'$\ \mu m$')
        plt.plot(x*1e6, Gauss([p[0],p[2]],x), 'k--')
        plt.plot(x*1e6, Gauss([p[1],p[3]],x), 'c--')
        plt.xlabel(r'$x\mathrm{\ [\mu m]}$')
        plt.ylabel('Counts')
        plt.ylim(bottom = 0.1)
        plt.ylim(top = max(left[0])*6)
        plt.legend()
        
        plt.subplot(122)
        right = plt.hist(ptcls[:,0]*1e6, bins = 50, weights = weights, log=False)
        plt.plot(x*1e6, fxg, label="G+G: "+r'$\sigma_1=$'+("%0.2f"%(np.abs(p[0])*1e6))+r'$\ \mu m$'+" & "+r'$\sigma_2=$'+("%0.2f"%(np.abs(p[1])*1e6))+r'$\ \mu m$')
        plt.plot(x*1e6, Gauss([p[0],p[2]],x), 'k--')
        plt.plot(x*1e6, Gauss([p[1],p[3]],x), 'c--')
        plt.xlabel(r'$x\mathrm{\ [\mu m]}$')
        plt.ylim(bottom = 0.1)
        plt.ylim(top = max(right[0])*1.2)
        plt.legend(); plt.show()
        
        GaussPlusGauss_Percent(p)
        
        return

class GaussianElectronBeam(ElectronBeam):
    """ A electron beam with a Gaussian transverse profile. 
    
    Parameters
    ----------
    gamma : double
        The relativistic factor of the beam.
    emittance : double
        The normalized emittance of the beam in m*rad.
    betax : double
        The beta function in the x-direction at z=0, in m.
    betay : double
        The beta function in the y direction at z=0, in m.
    alphax : double
        The alpha parameter of the beam in the x-direction at z=0.
    alphay : double
        The alpha parameter of the beam in the y-direction at z=0.
    sigmaz : double
        The RMS size of the beam in z in m.
    dE : double
        The energy spread of the beam as a fraction, 0.01 = +-1% energy spread.
    """
    
    def __init__(self, params):
        self.keys = self.keys.copy()
        self.keys.extend(
                ['gamma',
                 'emittance',
                 'betax',
                 'betay',
                 'alphax',
                 'alphay',
                 'sigmaz',
                 'dE'])
        super().__init__(params)
        
    def action_angle_distribution(self):
        """ Initialize particles in action-angle coordinates. 
        
        Returns
        -------
        ux : array of double
            Particle positions in ux.
        vx : array of double
            Particle positions in vx.
        uy : array of double
            Particle positions in uy.
        vy : array of double
            Particle positions in vy.
        """
        N = self.N
        gamma = self.gamma
        emittance = self.emittance
        # Calculate arrays of random numbers
        x1r = np.random.uniform(0, 1, N)
        x2r = np.random.uniform(0, 1, N)
        y1r = np.random.uniform(0, 1, N)
        y2r = np.random.uniform(0, 1, N)
        # Choose the particles from a distribution
        Jx = -emittance * np.log(x1r) / gamma
        Jy = -emittance * np.log(y1r) / gamma
        phix = 2*np.pi*x2r
        phiy = 2*np.pi*y2r
        ux = np.sqrt(2*Jx)*np.cos(phix)
        vx = -np.sqrt(2*Jx)*np.sin(phix)
        uy = np.sqrt(2*Jy)*np.cos(phiy)
        vy = -np.sqrt(2*Jy)*np.sin(phiy)
        return ux, vx, uy, vy
    
    def initialize_particles(self, offset_x=0.0, offset_y=0.0, offset_xp=0.0, offset_yp=0.0):
        """ Initialize the particles in a 6D distribution. """
        N = self.N
        gamma = self.gamma
        betax = self.betax
        betay = self.betay
        ptcls = np.zeros((N, 6), dtype='double')
        ux, vx, uy, vy = self.action_angle_distribution()
        # Calculate the coordinates
        ptcls[:, 0] = ux*np.sqrt(betax) + offset_x
        ptcls[:, 1] = (vx-self.alphax*ux) / np.sqrt(betax) + offset_xp
        ptcls[:, 2] = uy*np.sqrt(betay) + offset_y
        ptcls[:, 3] = (vy-self.alphay*uy) / np.sqrt(betay) + offset_yp
        ptcls[:, 4] = np.random.normal(0.0, self.sigmaz, N)
        ptcls[:, 5] = gamma * (1 + self.dE*np.random.uniform(-1, 1, N))
        super().initialize_particles(ptcls) 

class OffsetGaussianElectronBeam(GaussianElectronBeam):
    def __init__(self, params):
        self.keys = self.keys.copy()
        self.keys.extend(
                ['offset_x',
                 'offset_y',
                 'offset_xp',
                 'offset_yp'])
        super().__init__(params)
        
    def initialize_particles(self):
        """ Initialize the particles in a 6D distribution. """
        super().initialize_particles(self.offset_x, self.offset_y, self.offset_xp, self.offset_yp) 

class VorpalElectronBeam(ElectronBeam):
    """ A electron beam imported from VSim, includes weights. 
    
    Parameters
    ----------
    filename : string
        Filename of h5 file to load data from
    thresh : double
        Lowest particle weight to load, set to 0 to load all
    """
    
    def __init__(self, params):
        self.keys = self.keys.copy()
        self.keys.extend([
                'filename',
                'thresh',
                'minz',
                'maxz'])
        super().__init__(params)
    
    def initialize_particles(self):
        """ Initialize the particles in a 6D distribution. """
        #Given a filename, parse the info
        # 'Drive_Witness_Ramps_WitnessBeam_2.h5'
        filename = self.filename
        thresh = self.thresh
        minz = self.minz
        maxz = self.maxz
        
        file = filename.split("/")[-1]
        strlist = file.split("_")
        dumpInd = int(strlist[-1].split(".")[0])
        species = strlist[-2]
        simname = "_".join(strlist[:-2])
        
        data = Vload.get_species_data(filename, species)
        dim = int(data.attrs['numSpatialDims'])
        
        y = np.array(Vanalyze.get_y(data,dim)) #tran
        uy = np.array(Vanalyze.get_uy(data,dim))
        z = np.array(Vanalyze.get_z(data,dim)) #tran
        uz = np.array(Vanalyze.get_uz(data,dim))
        x = np.array(Vanalyze.get_x(data,dim)) #long
        gamma = np.array(Vanalyze.get_ptc_gamma(data))
        weights = np.array(Vanalyze.get_weights(data))
        
        ux = np.array(Vanalyze.get_ux(data,dim))
        uy = uy/ux
        uz = uz/ux
        
        x = x-np.average(x, weights=weights)#recenter to x=-0
        
        if minz > -np.inf:
            zset = np.array(np.where(x > minz)[0])
            y = y[zset]
            uy = uy[zset]
            z = z[zset]
            uz = uz[zset]
            x = x[zset]
            gamma = gamma[zset]
            weights = weights[zset]
            
        if maxz < np.inf:
            zset = np.array(np.where(x < maxz)[0])
            y = y[zset]
            uy = uy[zset]
            z = z[zset]
            uz = uz[zset]
            x = x[zset]
            gamma = gamma[zset]
            weights = weights[zset]
        
        threshset = np.array(np.where(weights > thresh)[0])
        N = len(threshset)
        self.N = N
        
        ptcls = np.zeros((N, 7), dtype='double')
        ptcls[:, 0] = y[threshset]
        ptcls[:, 1] = uy[threshset]
        ptcls[:, 2] = z[threshset]
        ptcls[:, 3] = uz[threshset]
        ptcls[:, 4] = x[threshset]
        ptcls[:, 5] = gamma[threshset]
        ptcls[:, 6] = weights[threshset]
        
        super().initialize_particles(ptcls)
    
    def get_weights(self, ptcls):
        return ptcls[:, 6]

    def get_emittance(self, ind):
        ptcls = self.load_ptcls(ind)[0]
        weights = self.get_weights(ptcls)
        return super().get_emittance(ind, weights)
        
    def get_gamma_n(self, ind):
        ptcls = self.load_ptcls(ind)[0]
        weights = self.get_weights(ptcls)
        return super().get_gamma_n(ind, weights)

    def get_emittance_n(self, ind):
        """ Calculate the normalized emittance from a particular save file. """
        ptcls = self.load_ptcls(ind)[0]
        weights = self.get_weights(ptcls)
        ex, ey = super().get_emittance(ind, weights=weights)
        gamma = super().get_gamma_n(ind, weights=weights)
        ex = ex*gamma
        ey = ey*gamma
        return ex, ey

    def get_sigmar(self, ind):
        ptcls = self.load_ptcls(ind)[0]
        weights = self.get_weights(ptcls)
        return super().get_sigmar(ind, weights)

    def get_sigmarp(self, ind):
        ptcls = self.load_ptcls(ind)[0]
        weights = self.get_weights(ptcls)
        return super().get_sigmarp(ind, weights)
    
    def plot_current_phase(self, xlim=None, ylim=None):
        weights = self.get_weights(self.ptcls)
        super().plot_phase(self.ptcls, self.z[-1], xlim, ylim, weights=weights)
        plt.show()
    
    def plot_phase_at(self, ind):
        ptcls, z = self.load_ptcls(ind)
        weights = self.get_weights(ptcls)
        super().plot_phase(ptcls, z, weights=weights)
        plt.show()
    
    def plot_phase(self, ptcls, z, xlim=None, ylim=None):
        weights = self.get_weights(ptcls)
        super().plot_phase(ptcls, z, xlim, ylim, weights)
        
    def plot_hist_at(self, ind):
        ptcls, z = self.load_ptcls(ind)
        weights = self.get_weights(ptcls)
        super().plot_hist(ptcls, z, ind, weights=weights)
        plt.show()
    
    def plot_hist(self, ptcls, z, ind):
        weights = self.get_weights(ptcls)
        super().plot_hist(ptcls, z, weights)
    
    def get_emittance_n_zcond(self, ind, minz, maxz):
        """ Calculate the normalized emittance from a particular save file. """
        ptcls = self.load_ptcls(ind)[0]
        weights = self.get_weights(ptcls)
        ex, ey = self.get_emittance_zcond(ind, minz, maxz, weights=weights)
        gamma = super().get_gamma_n(ind, weights=weights)
        ex = ex*gamma
        ey = ey*gamma
        return ex, ey    
    
    def get_emittance_zcond(self, ind, minz, maxz, weights, ptcls=None):
        """ Calculate the emittance from a particular save file. """
        ptcls = self.load_ptcls(ind)[0]
        x = self.get_x(ptcls)
        xp = self.get_xp(ptcls)
        y = self.get_y(ptcls)
        yp = self.get_yp(ptcls)
        #If weights aren't given, just initialize an array of 1's
            
        z = self.get_z(ptcls)
        zset = np.array(np.where((z > minz))[0])
        x = x[zset]
        xp = xp[zset]
        y = y[zset]
        yp = yp[zset]
        weights = weights[zset]
        z = z[zset]
        
        zset = np.array(np.where((z < maxz))[0])
        x = x[zset]
        xp = xp[zset]
        y = y[zset]
        yp = yp[zset]
        weights = weights[zset]
        
        # Calculate the differences from the average
        dx = x - np.average(x, weights=weights)
        dxp = xp - np.average(xp, weights=weights)
        dy = y - np.average(y, weights=weights)
        dyp = yp - np.average(yp, weights=weights)
        # Calculate the RMS sizes and the correlation
        sigmax2 = np.average(dx**2, weights=weights)
        sigmaxp2 = np.average(dxp**2, weights=weights)
        sigmaxxp = np.average(dx*dxp, weights=weights)
        sigmay2 = np.average(dy**2, weights=weights)
        sigmayp2 = np.average(dyp**2, weights=weights)
        sigmayyp = np.average(dy*dyp, weights=weights)
        # Calculate the emittance
        ex = np.sqrt(sigmax2*sigmaxp2 - sigmaxxp**2)
        ey = np.sqrt(sigmay2*sigmayp2 - sigmayyp**2)
        return ex, ey
    
#p[sig,gam,n0]
def Voigt(p, x):
    lorentz = (p[1]/(np.square(x)+np.square(p[1])))/np.pi
    gauss = np.exp(-.5*np.square(x)/np.square(p[0]))/p[0]/np.sqrt(2*np.pi)
    return p[2]*lorentz*gauss

def GaussPlusLorentz(p, x):
    lorentz = (p[1]/(np.square(x)+np.square(p[1])))/np.pi
    gauss = np.exp(-.5*np.square(x)/np.square(p[0]))
    return p[2]*gauss+p[3]*lorentz

def GaussPlusGauss(p, x):
    gauss1 = np.exp(-.5*np.square(x)/np.square(p[0]))
    gauss2 = np.exp(-.5*np.square(x)/np.square(p[1]))
    return p[2]*gauss1 + p[3]*gauss2

def GaussTimesGauss(p, x):
    gauss1 = np.exp(-.5*np.square(x)/np.square(p[0]))
    gauss2 = np.exp(-.5*np.square(x)/np.square(p[1]))
    return p[2]*gauss1 * gauss2

#p[gam,B]
def Lorentz(p,x):
    return p[1]*(p[0]/(np.square(x)+np.square(p[0])))/np.pi

#p[sig,A]
def Gauss(p, x):
    return p[1]*np.exp(-.5*np.square(x)/np.square(p[0]))
    
def FitDataSomething(data, axis, function, guess = [0.,0.,0.]):
    errfunc = lambda p, x, y: function(p, x) - y
    p0 = guess
    p1, success = optimize.leastsq(errfunc,p0[:], args=(axis, data))
    """
    plt.plot(axis, data, label=datlabel)
    #plt.plot(axis, function(guess,axis), label="Guess "+ function.__name__ +" profile")
    plt.plot(axis, function(p1,axis), label="Fitted "+ function.__name__ +" profile")
    plt.title("Comparison of data with "+ function.__name__ +" profile")
    plt.legend(); plt.grid(); plt.show()
    """
    return p1

#Given parameters for the Gauss + Gauss fit, what percentage is particles are in each Gaussian
def GaussPlusGauss_Percent(p):
    pi = [np.abs(p[0]), np.abs(p[2])]
    po = [np.abs(p[1]), np.abs(p[3])]
    inner = Int.quad(lambda x: Gauss(pi, x), -5*pi[0], 5*pi[0])[0]
    outer = Int.quad(lambda x: Gauss(po, x), -5*po[0], 5*po[0])[0]
    print(inner, outer)
    print("Inner: ",inner/(inner+outer)*100,"%")
    print("Outer: ",outer/(inner+outer)*100,"%")
    return    
    
    
    
    
    
    
    
    
    
    
    
    