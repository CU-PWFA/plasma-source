#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 16:57:52 2017

@author: litos
"""

def calc_Bmag(Tb,Tl):
    """ Calculates Bmag.
    
        Bmag measures matching of a beam to a lattice.
    
        Parameters
        ----------
        Tb
            Twiss parameter array for beam:
                Tb[0] = beta (m)
                Tb[1] = alpha
                Tb[2] = gamma (1/m)
        Tl
            Twiss parameters for lattice:
                Tl[0] = beta (m)
                Tl[1] = alpha
                Tl[2] = gamma (1/m)

        Returns
        -------
        Bmag
            matching parameter >= 1
                =1 -> matched beam
                >1 -> mismatched beam
    """
    [beta_b,alpha_b,gamma_b] = Tb
    [beta_l,alpha_l,gamma_l] = Tl
    
    Bmag = (1/2)*(beta_b*gamma_l-2*alpha_b*alpha_l+gamma_b*beta_l)

    return Bmag
