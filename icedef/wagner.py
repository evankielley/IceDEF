"""Iceberg drift model from Wagner et al, 2017.
"""

import numpy as np
import cmath
import netCDF4 as nc


def drift(iceberg, vau, vav, vwu, vwv, dt):
    
    # Constants
    R = 6378*1e3
    om = 7.2921e-5
    rhow = 1027
    rhoa = 1.2
    rhoi = iceberg.rho
    drho = rhow - rhoi
    Cw = iceberg.Cdw
    Ca = iceberg.Cda
    gam = np.sqrt(rhoa*drho/rhow/rhoi*(Ca/Cw))
   
    # Iceberg attributes
    t = iceberg.T  # time of the iceberg (datetime)
    x = iceberg.X  # x-component of iceberg position (degrees longitude)
    y = iceberg.Y  # y-component of iceberg position (degrees latitiude)
    l = iceberg.L  # length of the iceberg (m)
    w = iceberg.W  # width of the iceberg (m)
    
    
    # Drift
    S = np.pi*((l*w)/(l+w))
    ff = 2*om*np.sin((np.abs(y)*np.pi)/180)
    lam = np.sqrt(2)*Cw*(gam*np.sqrt(vau**2 + vav**2))/(ff*S)

    
    if lam < 0.1:
        alpha = lam*(lam**4*(lam**4*(lam**4*(-0.0386699020961393*lam**4 + \
            0.055242717280199) - 0.0883883476483184) + \
            0.176776695296637) - 0.707106781186548)
    else:
        alpha = np.multiply(np.divide(np.sqrt(2),np.power(lam, 3)),(1-np.sqrt(1+np.power(lam,4))))
        
    if lam < 0.6:
        beta = lam**3*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*\
            (0.0153268598203613*lam**4 - 0.0151656272365985) + \
            0.0180267866272764) + 0.0219176256311202) - \
            0.0274446790511418) + 0.0357675015202851) - \
            0.0493731785691779) + 0.0745776683282687) - \
            0.132582521472478) + 0.353553390593274)
    else:
        beta = np.real(np.multiply(np.divide(1.,np.power(lam,3.)),cmath.sqrt(np.multiply((4.+np.power(lam,4.)), \
            cmath.sqrt(1.+np.power(lam,4.)))-3.*np.power(lam,4.)-4.)))

    # Iceberg velocity   
    viu = vwu + gam*(-alpha*vav + beta*vau)
    viv = vwv + gam*(alpha*vau + beta*vav)
    
    return viu, viv

    
def melt(iceberg, vau, vav, vwu, vwv, sst, dt):
    
    # Constants
    sst0 = -4
    Cs1 = 1.5; Cs2 = 0.5; Cs3 = 0.1
    CMv1 = 7.62e-3; CMv2 = 1.29e-3; CMe1 = 0.5
    CMb1 = 0.58; CMb2 = 0.8; CMb3 = 0.2
    
    # Iceberg attributes
    t = iceberg.T  # time of the iceberg (datetime)
    x = iceberg.X  # x-component of iceberg position (degrees longitude)
    y = iceberg.Y  # y-component of iceberg position (degrees latitiude)
    l = iceberg.L  # length of the iceberg (m)
    w = iceberg.W  # width of the iceberg (m)
    h = iceberg.H  # height of the iceberg (m)
    

    # Melt Rates
    Me = CMe1*(Cs1*np.sqrt(vau**2 + vav**2)**Cs2 + Cs3*np.sqrt(vau**2 + vav**2))
    Mv = CMv1*sst + CMv2*sst**2
    Mb = CMb1*np.power(np.sqrt(np.square(viu-vwu)+np.square(viv-vwv)),CMb2)*(sst - sst0)/l**CMb3

    # Iceberg dimensions after melting
    l_new = l - (Mv + Me)*(dt/(24*3600))  # convert dt from secs to days
    w_new = w - (Mv + Me)*(dt/(24*3600))
    h_new = h - Mb*(dt/(24*3600))

    if w_new < 0.85*h_new:
        # Rollover
        print('rollover')
        w_new, h_new = h_new, w_new

    if w_new > l_new:
        # Ensure l is greater than w
        print('swap l and w')
        w_new, l_new = l_new, w_new
        
    return l_new, w_new, h_new
