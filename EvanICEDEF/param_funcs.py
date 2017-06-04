import numpy as np

def a(U):
    # \alpha in the paper
    a = np.multiply(np.divide(np.sqrt(2),np.power(U, 3)),(1-np.sqrt(1+np.power(U,4))))
    return a

def b(U):
    # \beta in the paper
    #print(np.multiply((4+np.power(U,4)),np.sqrt(1+np.power(U,4)))-3*np.power(U,4)-4)
    b = np.real(np.multiply(np.divide(1,np.power(U,3)),np.sqrt(np.multiply((4+np.power(U,4)), \
                                            np.sqrt(1+np.power(U,4)))-3*np.power(U,4)-4)))
    if not b >= 0:
        b = 0
    return b


def S(l, w):
    # Harmonic mean length
    #S = np.pi*(np.multiply(l, np.divide(w,(l+w))))
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
