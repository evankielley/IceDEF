"""Small mathematical functions necessary for calculating iceberg drift."""
import math
import cmath
import numpy as np
import sympy as sp

def alpha(U):
    if U < 0.1:
        a = U*(U**4*(U**4*(U**4*(-0.0386699020961393*U**4 + \
            0.055242717280199) - 0.0883883476483184) + \
            0.176776695296637) - 0.707106781186548)
    else:
        a = np.multiply(np.divide(np.sqrt(2),np.power(U, 3)),(1-np.sqrt(1+np.power(U,4))))
    return a
        
def beta(U):
    if U < 0.6:
        b = U**3*(U**4*(U**4*(U**4*(U**4*(U**4*(U**4*(U**4*(U**4*\
            (0.0153268598203613*U**4 - 0.0151656272365985) + \
            0.0180267866272764) + 0.0219176256311202) - \
            0.0274446790511418) + 0.0357675015202851) - \
            0.0493731785691779) + 0.0745776683282687) - \
            0.132582521472478) + 0.353553390593274)
    else:
        b = np.real(np.multiply(np.divide(1.,np.power(U,3.)),cmath.sqrt(np.multiply((4.+np.power(U,4.)), \
            cmath.sqrt(1.+np.power(U,4.)))-3.*np.power(U,4.)-4.)))
    return b
                
def S(l, w):
    """ Computes harmonic mean length"""
    S = np.pi*(np.divide(np.multiply(l,w),l+w))
    return S

def ff(lati,om):
    """ Computes latitude in degrees"""
    ff = 2*om*np.sin(np.abs(lati)*np.pi/180)
    ff = np.float64(ff)
    return ff

def Ut(u, lati, S, Cw, g, om):
    """Computes Lambda in the papers"""
    u = np.float64(u)
    lati = np.float64(lati)
    S = np.float64(S)
    Cw = np.float64(Cw)
    g = np.float64(g)
    om = np.float64(om)
    Ut = np.sqrt(2)*Cw*np.divide(g,ff(lati,om))*np.divide(u,S)
    return Ut

def find_nearest(array,value):
    """Returns the indice of the closest Lat or Lon to input y or x"""
    value = float(value)
    idx = (abs(array-value)).argmin()
    return idx

def find_berg_grid(LAT,LON,x,y):
    """Finds the tightest bracket of lat and lon for inputs y and x"""
    xi2a = np.where(LON <= x)[0][-1]
    xi2b = np.where(LON > x)[0][0]
    yi2a = np.where(LAT <= y)[0][-1]
    yi2b = np.where(LAT > y)[0][0]
    return xi2a,xi2b,yi2a,yi2b
