import math
import cmath
import numpy as np

def a(U):
    # \alpha in the paper
    a = np.multiply(np.divide(np.sqrt(2),np.power(U, 3)),(1-np.sqrt(1+np.power(U,4))))
    return a

def b(U):
    b = float(np.real(np.multiply(np.divide(1,np.power(U,3)),cmath.sqrt(np.multiply((4+np.power(U,4)), \
                                           cmath.sqrt(1+np.power(U,4)))-3*np.power(U,4)-4))))
    if math.isnan(b):
        b = 0
    return b

def S(l, w):
    # Harmonic mean length
    S = np.pi*(np.divide(np.multiply(l,w),l+w))
    return S

def ff(lati,om):
    # Latitude in degrees
    ff = 2*om*np.sin(abs(lati)*np.pi/180)
    return ff

def Ut(u, lati, S, Cw, g, om):
    # \Lambda in the papers
    Ut = np.sqrt(2)*Cw*g/ff(lati,om)*u/S
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
