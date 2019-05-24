#!/usr/bin/python

#============================
# make plasma density profile
#============================

## common libraries
from numpy import *

## define function
def makeplasma(np0,shape,z,args):

    if shape.lower() == 'gauss': # normal gauss. func.
        [z1,sig1,z2,sig2] = args
        npl = np0*((z<z1)*exp(-pow(z-z1,2)/(2*pow(sig1,2)))+\
                  (z>=z1)*(z<=z2)+\
                  (z>z2)*exp(-pow(z-z2,2)/(2*pow(sig2,2))))
    elif shape.lower() == 'gen_gauss': # generalized gauss. func.
        [z1,sig1,P1,sig2] = args
        z1   = float(z1)
        z2   = float(2*z[len(z)-1])
        sig1 = float(sig1)
        P1   = float(P1)
#        print('type P1: ',type(P1))
#        print('type z-z1: ',type(z-z1))
        
        npl = np0*((z<z1)*exp(-((z-z1)**P1)/(2*(sig1**P1)))+\
                   (z>=z1))
#                   (z>=z1)*(z<=z2)+\
#                   (z>z2)*exp(-((z-z2)**2)/(2*(sig2**2))))
        
        
#        npl = np0*((z<z1)*exp(-pow(z-z1,P1)/(2*pow(sig1,P1)))+\
#                   (z>=z1)*(z<=z2)+\
#                   (z>z2)*exp(-pow(z-z2,2)/(2*pow(sig2,2))))
    elif shape.lower() == 'trap': # trapezoidal func.
        [z1,sig1,z2,sig2] = args
        slope = (np0/2)/(z1-sig1)
        npl = np0*((z<=z1-2*sig1)*0+\
                   (z>z1-2*sig1)*(z<z1)*(1/(2*sig1))*(z-(z1-2*sig1))+\
                   (z>=z1)*(z<=z2)+\
                   (z>z2)*(-1/(2*sig1))*(z-(z1-2*sig1)))
    elif shape.lower() == 'sigmoid': # sigmoidal func.
        [z1,sig1,z2,sig2] = args
        npl = np0*(1/(1+exp(-(z-(z1))/sig1)))
    elif shape.lower() == 'gomp': # gompertz func.
        [z1,sig1,z2,sig2] = args
        npl = np0*(1-exp(-exp((z-z1)/sig1)))
        
    else:
        print('bad plasma ramp shape' ) 

    return npl
