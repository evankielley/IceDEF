"""Iceberg drift model from Wagner et al, 2017.

This module can simulate the drifting and melting of an iceberg over one timestep.

For information on the equations of motion used, please see:
    Wagner, T.J., R.W. Dell, and I. Eisenman, 2017: An Analytical Model of Iceberg Drift.
    J. Phys. Oceanogr., 47, 1605–1616, https://doi.org/10.1175/JPO-D-16-0262.1
"""

import numpy as np
import cmath


def drift(t, x, y, vx, vy, C):
    """This function simulates the drift of an iceberg over one timestep."""

    om = C['OM']
    rhoa = C['RHOA']
    rhow = C['RHOW']
    rhoi = C['RHOI']
    Cw = C['Cdw']
    Ca = C['Cda']
    l = C['l']
    w = C['w']
    vax = C['vax']
    vay = C['vay']
    vwx = C['vwx']
    vwx = C['vwy']

    drho = rhow - rhoi
    gam = np.sqrt((rhoa*drho)/(rhow*rhoi)*(Ca/Cw))  # dimensionless parameter


    # Drift
    S = np.pi*((l*w)/(l+w))
    ff = 2*om*np.sin((np.abs(y)*np.pi)/180)
    lam = np.sqrt(2)*Cw*(gam*np.sqrt(vax**2 + vay**2))/(ff*S)


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
    vx = vwx + gam*(-alpha*vay + beta*vax)
    vy = vwy + gam*(alpha*vax + beta*vay)

    return vx, vy


def melt(dt, l, w, h, vx, vy, C):
    """This function simulates the melting of an iceberg after one timestep."""

    sst0 = C['sst0']
    Cs1 = C['Cs1']
    Cs2 = C['Cs2']
    Cs3 = C['Cs3']
    CMv1 = C['CMv1']
    CMv2 = C['CMv2']
    CMe1 = C['CMe1']
    CMb1 = C['CMb1']
    CMb2 = C['CMb2']
    CMb3 = C['CMb3']
    vax = C['vax']
    vay = C['vay']
    vwx = C['vwx']
    vwx = C['vwy']
    sst = C['sst']

    # Example C
    #sst0 = -4
    #Cs1 = 1.5; Cs2 = 0.5; Cs3 = 0.1
    #CMv1 = 7.62e-3; CMv2 = 1.29e-3; CMe1 = 0.5
    #CMb1 = 0.58; CMb2 = 0.8; CMb3 = 0.2


    # Melt Rates
    Me = CMe1*(Cs1*np.sqrt(vax**2 + vay**2)**Cs2 + Cs3*np.sqrt(vax**2 + vay**2))
    Mv = CMv1*sst + CMv2*sst**2
    Mb = CMb1*np.power(np.sqrt(np.square(vx-vwx)+np.square(vy-vwy)),CMb2)*(sst - sst0)/l**CMb3

    # Iceberg dimensions after melting
    l -= (Mv + Me)*(dt/(24*3600))  # convert dt from secs to days
    w -= (Mv + Me)*(dt/(24*3600))
    h -= Mb*(dt/(24*3600))

    if w < 0.85*h:
        # Rollover
        print('rollover')
        w, h = h, w

    if w > l:
        # Ensure l is greater than w
        print('swap l and w')
        w, l = l, w

    return l, w, h
