# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 11:56:55 2022

@author: Robert
"""

import numpy as np
from PIL import Image
from scipy import ndimage
import base64
import ast
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from scipy.optimize import curve_fit
import libtiff
import glob
import os
import re

PATH = ''

# TODO - rewrite this as a class like the image stuff is
def load_PHASE(dataset, name):
    dataset = str(dataset)
    year, month, day = get_date_from_dataset(dataset)
    filename = PATH + 'WAVEFRONT' + '/year_{}/month_{}/day_{}'.format(year, month, day) + '/' + name + '.txt'
    data = np.loadtxt(filename, dtype='S', delimiter='\t', skiprows=8)
    sel = data == b''
    data[sel] = 'Nan'
    return np.array(data, dtype='float')

def plot_phase(image):
    im_max = np.nanmax(image)
    im_min = np.nanmin(image)
    cen = 0.5*(im_max-im_min)+im_min
    fig = plt.figure(figsize=(4.85, 3), dpi=300)
    ax = plt.subplot()
    im = ax.imshow(image-cen, cmap='RdYlBu')#, vmin=-im_max, vmax=im_max)
    cb = fig.colorbar(im)
    cb.set_label(r'Phase ($\mathrm{\mu m}$)')
    ax.set_xlabel(r'$x$ (px)')
    ax.set_ylabel(r'$y$ (px)')
    ax.set_xlim(-0.5, image.shape[1]-0.5)
    ax.set_ylim(image.shape[0]-0.5, -0.5)
    plt.text(0, 2, r"P-V Phase: %0.3f $\mathrm{\mu m}$"%(im_max-im_min))
    plt.text(0, 4, r"RMS Phase: %0.3f $\mathrm{\mu m}$"%(np.std(image, where=~np.isnan(image))))
    return fig, ax, im, cb

def plot_phase_error(image, poly, subtracted=False):
    tilt_h = poly.polynomial_phase_haso(1)
    tilt_v = poly.polynomial_phase_haso(2)
    focus = poly.polynomial_phase_haso(3)
    if subtracted:
        fig, ax, im, cb = plot_phase(image)
    else:
        fig, ax, im, cb = plot_phase(image-tilt_h-tilt_v-focus)
    plt.text(0, 6, r"Tilt H: %0.3f $\mathrm{\mu m}$"%(poly.coefficients[0]))
    plt.text(0, 8, r"Tilt V: %0.3f $\mathrm{\mu m}$"%(poly.coefficients[1]))
    plt.text(0, 10, r"Focus : %0.3f $\mathrm{\mu m}$"%(poly.coefficients[2]))
    return fig, ax, im, cb

class Zernike():
    rexp = {
        'sensor' : re.compile(r'ASO Size : (?P<X>[0-9]{2}) x (?P<Y>[0-9]{2}); ASO Step : (?P<dx>.*) um x (?P<dy>.*) um'),
        'pupil' : re.compile(r'Pupil Center \(x, y\) : \((?P<pX>.*) mm x (?P<pY>.*) mm\) ; Pupil radius : (?P<pR>.*) mm;')
        }
    
    def __init__(self, dataset, name):
        self.dataset = str(dataset)
        year, month, day = get_date_from_dataset(dataset)
        self.filename = get_path('WAVEFRONT', year, month, day) + name + '.txt'
        data = np.loadtxt(self.filename, dtype='S', delimiter='\t', skiprows=6)
        # The P-V coefficients are simply dump multiplication factors on the polynomials
        self.coefficients = np.zeros(32)
        for i in range(32):
            self.coefficients[i] = float(data[i, 3])
        with open(self.filename, 'r') as f:
            for line in f:
                key, match = parse_line(line, self.rexp)
                if key == 'sensor':
                    self.X = int(match.group('X'))
                    self.Y = int(match.group('Y'))
                    self.dx = float(match.group('dx'))
                    self.dy = float(match.group('dy'))
                if key == 'pupil':
                    self.pX = float(match.group('pX'))
                    self.pY = float(match.group('pY'))
                    self.pR = float(match.group('pR'))
        self.poly_enum = {
            1: self.tilt0,
            2: self.tilt90,
            3: self.focus,
            4: self.astigmatism0,
            5: self.astigmatism45,
            6: self.coma0,
            7: self.coma90,
            8: self.spherical3rd,
            9: self.trefoil0,
            10: self.trefoil90,
            11: self.astigmatism5th0,
            12: self.astigmatism5th45,
            13: self.coma5th0,
            14: self.coma5th90,
            15: self.spherical5th,
            16: self.tetrafoil0,
            17: self.tetrafoil45,
            18: self.trefoil7th0,
            19: self.trefoil7th90,
            20: self.astigmatism7th0,
            21: self.astigmatism7th45,
            22: self.coma7th0,
            23: self.coma7th90,
            24: self.spherical7th,
            25: self.pentafoil0,
            26: self.pentafoil90,
            27: self.tetrafoil9th0,
            28: self.tetrafoil9th45,
            29: self.trefoil9th0,
            30: self.trefoil9th90,
            31: self.astigmatism9th0,
            32: self.astigmatism9th45
        }
    
    def polynomial_phase(self, ind, x, y):
        """ Calculate the phase of a given polynomial on the passed grid. x and y should be in mm. """
        r = np.sqrt(x[None, :]**2+y[:, None]**2)
        phi = np.arctan(y[:, None]/x[None, :])
        return self.poly_enum[ind](r, phi)
    
    def polynomial_phase_haso(self, ind):
        x = np.arange(self.X)*self.dx*1e-3
        y = np.arange(self.Y)*self.dy*1e-3
        x -= 0.5*x[-1]+self.pX
        y -= 0.5*y[-1]-self.pY
        r = np.sqrt(x[None, :]**2+y[:, None]**2)
        phi = -np.arctan(y[:, None]/x[None, :])
        sel = x<0
        phi[:, sel] -= np.pi
        sel = r > self.pR
        r[sel] = np.NAN
        return self.poly_enum[ind](r, phi)
    
    def build_phase(self, inds, x, y):
        """ Calculate the sum of polynomials on the passed grid. x and y should be in mm. """
        phase = np.zeros((len(y), len(x)))
        r = np.sqrt(x[None, :]**2+y[:, None]**2)
        phi = np.arctan(y[:, None]/x[None, :])
        for ind in inds:
            phase += self.poly_enum[ind](r, phi)
        return phase
    
    def tilt0(self, r, phi):
        return self.coefficients[0]*(r/self.pR)*np.cos(phi)
    
    def tilt90(self, r, phi):
        return self.coefficients[1]*(r/self.pR)*np.sin(phi)
    
    def focus(self, r, phi):
        return self.coefficients[2]*(2*(r/self.pR)**2-1)
    
    def astigmatism0(self, r, phi):
        return self.coefficients[3]*(r/self.pR)**2*np.cos(2*phi)
    
    def astigmatism45(self, r, phi):
        return self.coefficients[4]*(r/self.pR)**2*np.sin(2*phi)
    
    def coma0(self, r, phi):
        return self.coefficients[5]*(3*(r/self.pR)**2-2)*(r/self.pR)*np.cos(phi)
    
    def coma90(self, r, phi):
        return self.coefficients[6]*(3*(r/self.pR)**2-2)*(r/self.pR)*np.sin(phi)
    
    def spherical3rd(self, r, phi):
        return self.coefficients[7]*(6*(r/self.pR)**4-6*(r/self.pR)**2+1)
    
    def trefoil0(self, r, phi):
        return self.coefficients[8]*(r/self.pR)**3*np.cos(3*phi)
    
    def trefoil90(self, r, phi):
        return self.coefficients[9]*(r/self.pR)**3*np.sin(3*phi)
    
    def astigmatism5th0(self, r, phi):
        return self.coefficients[10]*(4*(r/self.pR)**2-3)*(r/self.pR)**2*np.cos(2*phi)
    
    def astigmatism5th45(self, r, phi):
        return self.coefficients[11]*(4*(r/self.pR)**2-3)*(r/self.pR)**2*np.sin(2*phi)
    
    def coma5th0(self, r, phi):
        return self.coefficients[12]*(10*(r/self.pR)**4-12*(r/self.pR)**2+3)*(r/self.pR)*np.cos(phi)
    
    def coma5th90(self, r, phi):
        return self.coefficients[13]*(10*(r/self.pR)**4-12*(r/self.pR)**2+3)*(r/self.pR)*np.sin(phi)
    
    def spherical5th(self, r, phi):
        return self.coefficients[14]*(20*(r/self.pR)**6-30*(r/self.pR)**4+12*(r/self.pR)**2-1)
    
    def tetrafoil0(self, r, phi):
        return self.coefficients[15]*(r/self.pR)**4*np.cos(4*phi)
    
    def tetrafoil45(self, r, phi):
        return self.coefficients[16]*(r/self.pR)**4*np.sin(4*phi)
    
    def trefoil7th0(self, r, phi):
        return self.coefficients[17]*(5*(r/self.pR)**2-4)*(r/self.pR)**3*np.cos(3*phi)
    
    def trefoil7th90(self, r, phi):
        return self.coefficients[18]*(5*(r/self.pR)**2-4)*(r/self.pR)**3*np.sin(3*phi)
    
    def astigmatism7th0(self, r, phi):
        return self.coefficients[19]*(15*(r/self.pR)**4-20*(r/self.pR)**2+6)*(r/self.pR)**2*np.cos(2*phi)
    
    def astigmatism7th45(self, r, phi):
        return self.coefficients[20]*(15*(r/self.pR)**4-20*(r/self.pR)**2+6)*(r/self.pR)**2*np.sin(2*phi)
    
    def coma7th0(self, r, phi):
        return self.coefficients[21]*(35*(r/self.pR)**6-60*(r/self.pR)**4+30*(r/self.pR)**2-4)*(r/self.pR)*np.cos(phi)
    
    def coma7th90(self, r, phi):
        return self.coefficients[22]*(35*(r/self.pR)**6-60*(r/self.pR)**4+30*(r/self.pR)**2-4)*(r/self.pR)*np.sin(phi)
    
    def spherical7th(self, r, phi):
        return self.coefficients[23]*(70*(r/self.pR)**8-140*(r/self.pR)**6+90*(r/self.pR)**4-20*(r/self.pR)**2+1)
    
    def pentafoil0(self, r, phi):
        return self.coefficients[24]*(r/self.pR)**5*np.cos(5*phi)
    
    def pentafoil90(self, r, phi):
        return self.coefficients[25]*(r/self.pR)**5*np.cos(5*phi)
    
    def tetrafoil9th0(self, r, phi):
        return self.coefficients[26]*(6*(r/self.pR)**2-5)*(r/self.pR)**4*np.cos(4*phi)
    
    def tetrafoil9th45(self, r, phi):
        return self.coefficients[27]*(6*(r/self.pR)**2-5)*(r/self.pR)**4*np.sin(4*phi)
    
    def trefoil9th0(self, r, phi):
        return self.coefficients[28]*(21*(r/self.pR)**4-30*(r/self.pR)**2+10)*(r/self.pR)**3*np.cos(3*phi)
    
    def trefoil9th90(self, r, phi):
        return self.coefficients[29]*(21*(r/self.pR)**4-30*(r/self.pR)**2+10)*(r/self.pR)**3*np.sin(3*phi)
    
    def astigmatism9th0(self, r, phi):
        return self.coefficients[30]*(56*(r/self.pR)**6-105*(r/self.pR)**4+60*(r/self.pR)**2-10)*(r/self.pR)**2*np.cos(2*phi)
    
    def astigmatism9th45(self, r, phi):
        return self.coefficients[31]*(56*(r/self.pR)**6-105*(r/self.pR)**4+60*(r/self.pR)**2-10)*(r/self.pR)**2*np.sin(2*phi)
    
