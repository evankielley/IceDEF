"""Iceberg drift model from Wagner, 2017.
.. module:: wagner
    :platform: Unix, Windows
    :synopsis: This module computes the trajectory of an iceberg in terms of it's longitude and latitude. The equations used in this
    model stem from the paper, An Analytical Model of Iceberg Drift, Wagner et. al. (2017).
.. moduleauthor:: Evan Kielley <evankielley@gmail.com>
"""

import numpy as np
import cmath

# Constants
R = 6378*1e3
om = 7.2921e-5
rhow = 1027
rhoa = 1.2
rhoi = 850
drho = rhow - rhoi
Cw = 0.9
Ca = 1.3
gam = np.sqrt(rhoa*drho/rhow/rhoi*(Ca/Cw))
sst0 = -4
Cs1 = 1.5; Cs2 = 0.5; Cs3 = 0.1
CMv1 = 7.62e-3; CMv2 = 1.29e-3; CMe1 = 0.5
CMb1 = 0.58; CMb2 = 0.8; CMb3 = 0.2



def wagner_drift(x, y, l, w, h, UA, VA, UW, VW, SST, t_ocean, t_atm, dt):
    """This functions computes the trajectory of an iceberg over one timestep.
    :param x: The iceberg's longitudinal coordinate (degrees).
    :type x: float.
    :param y: The iceberg's latitudinal coordinate (degrees).
    :type y: float.
    :param l: The iceberg's length dimension (m).
    :type l: float.
    :param w: The iceberg's width dimension (m).
    :type w: float.
    :param h: The iceberg's height dimension (m).
    :type h: float.
    :param UA: The interpolated u-component of the wind velocity (m/s), East is positive.
    :type UA: array_like, shape (t, x, y).
    :param VA: The interpolated v-component of the wind velocity (m/s), North is positive.
    :type VA: array_like, shape (t, x, y).
    :param UW: The interpolated u-component of the current velocity (m/s), East is positive.
    :type UW: array_like, shape (t, x, y).
    :param VW: The interpolated v-component of the current velocity (m/s), North is positive.
    :type VW: array_like, shape (t, x, y).
    :param SST: The interpolated sea-surface temperature (degrees C).
    :type SST: array_like, shape (t, x, y).
    :param t_ocean: The timestep that agrees with the time units for the ocean current field.
    :type t_ocean: float.
    :param t_ocean: The timestep that agrees with the time units for the wind field.
    :type t_atm: float.
    :param dt: The timestep length (s).
    :type dt: float.
    :return x_new: The new longitudinal position of the iceberg (degrees).
    :rtype x_new: float.
    :return y_new: The new latitudinal position of the iceberg (degrees).
    :rtype y_new: float.
    :return l_new: The new length of the iceberg (m).
    :rtype l_new: float.
    :return w_new: The new width of the iceberg (m).
    :rtype w_new: float.
    :return h_new: The new height of the iceberg (m).
    :rtype h_new: float.
    """


    # Extract values from input fields
    
    vau = UA([t_atm, y, x])[0]
    vav = VA([t_atm, y, x])[0]  
    vwu = UW([t_ocean, y, x])[0] 
    vwv = VW([t_ocean, y, x])[0]
    sst = SST([t_ocean, y, x])[0]
    
    
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

        
    viu = vwu + gam*(-alpha*vav + beta*vau)
    viv = vwv + gam*(alpha*vau + beta*vav)

    y_new = y + (viv*dt)*(180/(np.pi*R))
    x_new = x + (viu*dt)/(np.cos((((y + y_new)/2)*np.pi)/180))*(180/(np.pi*R))


    # Decay

    Me = CMe1*(Cs1*np.sqrt(vau**2 + vav**2)**Cs2 + Cs3*np.sqrt(vau**2 + vav**2))
    Mv = CMv1*sst + CMv2*sst**2
    Mb = CMb1*np.power(np.sqrt(np.square(viu-vwu)+np.square(viv-vwv)),CMb2)*(sst - sst0)/l**CMb3

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
        
    return x_new, y_new, l_new, w_new, h_new
