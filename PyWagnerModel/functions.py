import math
import cmath
import numpy as np
from decimal import Decimal

def a(U):
    # \alpha in the paper
    a = np.multiply(np.divide(np.sqrt(2),np.power(U, 3)),(1-np.sqrt(1+np.power(U,4))))
    return a

def b_small(U,mU,mBETA):
    b = np.real(np.multiply(np.power(U,3)/np.sqrt(8),cmath.sqrt(1 - 3/4*np.power(U,4) \
        + 9/16*np.power(U,8) - 7/16*np.power(U,12) + 45/256*np.power(U,16))))
    if math.isnan(b) or b==0 or b==0.0:
        print('python UT: {}, b: {}, matlab UT: {}, b: {}'.format(U,b,mU,mBETA))
        print('corrected beta')
        b = mBETA
    return b

def b_big(U,mU,mBETA):
    b = np.real(np.multiply(np.divide(1.,np.power(U,3.)),cmath.sqrt(np.multiply((4.+np.power(U,4.)), \
        cmath.sqrt(1.+np.power(U,4.)))-3.*np.power(U,4.)-4.)))
    if math.isnan(b) or b==0 or b==0.0:
        print('python UT: {}, b: {}, matlab UT: {}, b: {}'.format(U,b,mU,mBETA))
        print('corrected beta')
        b = mBETA
    return b

def S(l, w):
    # Harmonic mean length
    S = np.pi*(np.divide(np.multiply(l,w),l+w))
    return S

def ff(lati,om):
    # Latitude in degrees
    ff = 2*om*np.sin(np.abs(lati)*np.pi/180)
    ff = np.float64(ff)
    return ff

def Ut(u, lati, S, Cw, g, om):
    # \Lambda in the papers
    u = np.float64(u)
    lati = np.float64(lati)
    S = np.float64(S)
    Cw = np.float64(Cw)
    g = np.float64(g)
    om = np.float64(om)
    Ut = np.sqrt(2)*Cw*np.divide(g,ff(lati,om))*np.divide(u,S)
    return Ut

def find_nearest(array,value):
    value = float(value)
    idx = (abs(array-value)).argmin()
    return idx

def find_berg_grid(LAT,LON,x,y):
    xi2a = np.where(LON <= x)[0][-1]
    xi2b = np.where(LON > x)[0][0]
    yi2a = np.where(LAT <= y)[0][-1]
    yi2b = np.where(LAT > y)[0][0]
    return xi2a,xi2b,yi2a,yi2b
