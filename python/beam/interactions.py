#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 14:56:11 2017

@author: robert
"""

import numpy as np
import beam.beams as beams
import beam.elements as elements
import beam.calc.plasma as pcalc
import beam.calc.laser as lcalc
import beam.calc.electron as ecalc


def pulse_plasma(pulse, plasma):
    """ Propagates a pulse through a gas, ionizing and refracting as it goes.
    
    Parameters
    ----------
    pulse : Pulse class
        The laser pulse to propagate through the plasma.
    plasma : Plasma class
        The gas to propagate the laser pulse through.
    """
    pulse.e = pcalc.plasma_refraction(pulse.e, pulse.x, pulse.y,
                      plasma.z, pulse.t, pulse.lam, plasma.n0, pulse.z[-1],
                      pulse.fft, pulse.ifft, pulse.save_field, 
                      plasma.save_plasma_density, plasma.atom, 
                      plasma.load_num_den, plasma.load_plasma_den)


def beam_phase(beam, phase):
    """ Applies a phase mask to a optical beam, either a pulse or laser.
    
    Parameters
    ----------
    beam : Pulse or Laser class
        The optical beam the apply the phase mask to.
    phase : Phase class
        The phase mask to apply to the beam.
    """
    beam.set_field(beam.e * np.exp(1j*phase.phi))


def beam_plasma(beam, plasma):
    """ Propagates a weak beam through a plasma without an ionization.
    
    Parameters
    ----------
    beam : Laser class
        The optical beam to propagate through the plasma.
    plasma : Plasma class
        The plasma to propagate the laser pulse through.
    """
    nh = 1.0 + plasma.n0*plasma.atom['alpha']*5.0e-8
    def loadn(ind):
        ne = plasma.load_plasma_den(ind)
        nplasma = -beam.lam**2 * 4.47869e-5
        ngas = plasma.atom['alpha']*5.0e-8
        dn = nplasma - ngas
        return dn * ne
    beam.e = lcalc.beam_prop(beam.e, beam.x, beam.y, plasma.z, beam.lam, nh,
                             beam.z[-1], beam.fft, beam.ifft, beam.save_field,
                             loadn)


def electron_plasma(electron, plasma):
    """ Propagate an electron beam through an ion column. """
    # TODO add in energy gain/loss in the plasma
    electron.ptcls = ecalc.electron_propagation_plasma(electron.ptcls,
                            plasma.z, 0.0, plasma.ne,
                            electron.save_ptcls, plasma.dgammadz)
