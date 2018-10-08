import numpy as np
from icedef.constants import *


def newtonian_drift(iceberg_velocity, current_velocity, wind_velocity, **kwargs):
    """Computes instantaneous iceberg acceleration."""

    # Constants
    Omega = EARTH_ROTATION_RATE
    rhoa = AIR_DENSITY
    rhow = SEAWATER_DENSITY

    # Args
    Vx, Vy = iceberg_velocity
    Vwx, Vwy = wind_velocity
    Vcx, Vcy = current_velocity

    # Kwargs
    Amwx, Amwy = kwargs.pop('current_acceleration', (0, 0))
    Ca = kwargs.pop('form_drag_coefficient_in_air', 1.5)
    Cw = kwargs.pop('form_drag_coefficient_in_water', 1.5)
    Cda = kwargs.pop('skin_drag_coefficient_in_air', 2.5e-4)
    Cdw = kwargs.pop('skin_drag_coefficient_in_water', 5e-4)
    As = kwargs.pop('sail_area', 9600.0)
    Ak = kwargs.pop('keel_area', 48000.0)
    At = kwargs.pop('top_area', 25600.0)
    Ab = kwargs.pop('bottom_area', 25600.0)
    M = kwargs.pop('mass', 5468160000.0)
    phi = kwargs.pop('latitude', 50)

    # Wind force
    Fax = (0.5 * rhoa * Ca * As + rhoa * Cda * At) * np.sqrt((Vwx - Vx)**2 + (Vwy - Vy)**2) * (Vwx - Vx)
    Fay = (0.5 * rhoa * Ca * As + rhoa * Cda * At) * np.sqrt((Vwx - Vx)**2 + (Vwy - Vy)**2) * (Vwy - Vy)

    # Current force
    Fwx = (0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * np.sqrt((Vcx - Vx)**2 + (Vcy - Vy)**2) * (Vcx - Vx)
    Fwy = (0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * np.sqrt((Vcx - Vx)**2 + (Vcy - Vy)**2) * (Vcy - Vy)

    # Coriolis force
    f = 2 * Omega * np.sin(np.deg2rad(phi))
    Fcx = f * M * Vy
    Fcy = -f * M * Vx

    # Water pressure force
    Vmwx = Vcx
    Vmwy = Vcy
    Amwx = Amwx
    Amwy = Amwy
    Fwpx = M * (Amwx + f * Vmwy)
    Fwpy = M * (Amwy - f * Vmwx)

    # Iceberg acceleration
    Ax = (Fax + Fwx + Fcx + Fwpx) / M
    Ay = (Fay + Fwy + Fcy + Fwpy) / M
    
    forces = [Fax, Fay, Fwx, Fwy, Fcx, Fcy, Fwpx, Fwpy] 

    return Ax, Ay, forces
