'''
Module used for running a WARGSim simulation of a beam propagating through a 
PWFA

author : keenan
'''

# standard python code
import numpy as np
import sys
import scipy.constants as const
import psutil
me = const.physical_constants['electron mass energy equivalent in MeV'][0]
me = me * 1e6
# Add WARGSim directory to path
sys.path.insert(0, "/home/keenan/WARGSim/")
from beams import electronbeam
from beams import betatronbeam
import calc.electron as ecalc
from elements import pwfa as pwfa
from interactions import interactions



# Useful functions for checking emittance growth

def match_beta(np0,gb0):
    np0 = np0*(1e17)
    kbeta = (1.33e-4)*np.sqrt(np0/gb0)
    beta = 1.0/kbeta
    return beta

def calc_Bmag(beta,beta_m):
    return (1/2)*(beta/beta_m + beta_m/beta)

def calc_beta_fact(Bmag):
    return (Bmag+np.sqrt((Bmag**2) - 1))

def init_beam(beam_params, n0):
    '''
    Function to initialize electron beam
    
    Parameters:
    -----------
    beam_params : dictionary
        Dictionary of beam parameters containing:
        'N'       : number of particles
        'beamE'   : beam energy (eV)
        'eps_n0'  : normalized emittance (m-rad)
        'beta0'   : array of beta function at z=0 (x,y)
        'alpha0'  : array of alpha function at z=0 (x,y)
        'rms_z0'  : rms bunch length (m)
        'rms_gb0' : relative energy spread
        'path'    : path to store beam dumps
        'B_mag'   : target B_mag for the beam
    
    n0 : float
        PWFA flat-top density (normalized to 1e17/cc)
        
    Returns:
    --------
    my_ebeam : obj
        Instance of the ElectronBeam object
    '''
    
    # Get beam centroid gamma
    gb0 = beam_params['beamE'] / me
    # Calculate matched beta 
    beta_m = match_beta(n0, gb0)
    beta_x00 = calc_beta_fact(beam_params['B_mag'])*beta_m
    beta_y00 = calc_beta_fact(beam_params['B_mag'])*beta_m
    beta_x0  = beam_params['beta0'][0]
    Bmag = calc_Bmag(beta_x0,beta_m)
    print(f'beta factor = {beta_x0/beta_m :.2f}')
    print(f'Bmag = {Bmag :.2f}')
    
    electronParams = {
            'name'    : 'ElectronBeam',
            'path'    : beam_params['path'],
            'load'    : False,
            'N'       : beam_params['N'],
            'shape'   : 'Gauss',
            'gb0'     : gb0,
            'rms_gb'  : beam_params['rms_gb0'],
            'rms_z'   : beam_params['rms_z0'],
            'eps_nx'  : beam_params['eps_n0'],
            'eps_ny'  : beam_params['eps_n0'],
            'beta_x'  : beta_x00, 
            'beta_y'  : beta_y00,
            'alpha_x' : beam_params['alpha0'][0],
            'alpha_y' : beam_params['alpha0'][1]
            }

    my_ebeam = electronbeam.ElectronBeam(electronParams)
    eps_nx0, eps_ny0 = my_ebeam.get_emit_n()
    print('init eps_nx = ',eps_nx0)
    print('init eps_ny = ',eps_ny0)
    return my_ebeam

def init_plasma(beam_params, plasma_params, hw_dn = 0):
    '''
    Function to intialize the PWFA
    
    Parameters:
    -----------
    beam_params : dictionary
                  Dictionary of beam parameters defined above
    plasma_params : dictionary
        Dictionary of plasma parameters containing:
        n0    : Density of flat top (normalized to 1e17/cc)
        L_ft  : Length of flat top (m)
        hw_up : half-width of up ramp
        shape : shape of plasma 
    
    Returns:
    --------
    pwfa0 : obj
        Instance of PWFA class
    '''
    
    def match_hw(np0, gb0, beta0):
        kbeta = (1.33e-4)*np.sqrt(np0/gb0)
        beta = beta0*kbeta
        hw = (1.96e-3)*beta**2 + 1.49*beta - 3.08
        return hw/kbeta

    def match_zw(np0,gb0,beta0):
        kbeta = (1.33e-4)*np.sqrt(np0/gb0)
        beta = beta0*kbeta
        zw = (-1.47e-2)*beta**2 - 4.34*beta + 17.9
        return zw/kbeta
   
   # Calculate ramp lengths and full length of plasma
    L_up = 8 * plasma_params['hw_up']
    gb0  = beam_params['beamE'] / me
    if hw_dn != 0:
        hw_dn = 1.00 * match_hw(plasma_params['n0']*1e17, 2 * gb0,\
                           beam_params['beta0'][0])
    L_dn  = min(1.600, 8 * hw_dn)
    L_p   = L_up + plasma_params['L_ft'] + L_dn
    pwfaParams = {
        'name'   : 'PWFA0',
        'L'      : L_p,
        'n0'     : plasma_params['n0'],
        'autoNz' : True,
        'gb0'    : gb0,
        'shape'  : plasma_params['shape'],
        'L_ft'   : plasma_params['L_ft'],
        'L_up'   : plasma_params['L_up'],
        'L_dn'   : plasma_params['L_dn'],
        'hw_up'  : plasma_params['hw_up'],
        'hw_dn'  : plasma_params["hw_dn"]
}

    pwfa0 = pwfa.PWFA(pwfaParams)
    return pwfa0
    
def propagate(beam_params, plasma_params, dumpPeriod = 1e9, name = ""):
    '''
    Function to propagate an electron beam through a PWFA using WARGSim
    
    Parameters:
    -----------
    beam_params: dictionary
        Dictionary of beam parameters defined above
    plasma_params: dictionary
        Dictionary of plasma parameters defined above
    dumpPeriod : int, optional
        Dump period of the simulation
    name : str, optional
        Name of the simulation run
    '''
    # Assign dump path to default or default + name
    path = "/home/keenan/Documents/data/Wargsim/dumps/"
    path = path + name
    beam_params['path'] = path
    
    my_ebeam = init_beam(beam_params, plasma_params['n0'])
    pwfa0   = init_plasma(beam_params, plasma_params)
    print("Step size", pwfa0.dz)
    # Get number of threads
    Ncores = psutil.cpu_count()
    numThread  = Ncores # int, number of threads to run on
    
    ecalc.ebeam_prop_pwfa_dumps(my_ebeam.ptcls, pwfa0.dz, pwfa0.ne, pwfa0.dgb,
                                dumpPeriod, my_ebeam.save_ptcls, numThread)
