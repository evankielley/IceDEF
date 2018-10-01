import numpy as np
from icedef.constants import *


def newtonian_drift(iceberg_velocity, current_velocity, wind_velocity, iceberg_constants):
    """Computes instantaneous iceberg acceleration."""

    # Constants
    Omega = EARTH_ROTATION_RATE
    rhoa = AIR_DENSITY
    rhow = SEAWATER_DENSITY

    # Args
    Vx, Vy = iceberg_velocity

    Vwx, Vwy = wind_velocity
    Vcx, Vcy = current_velocity
    # Amwx, Amwy = current_acceleration
    Amwx, Amwy = 0, 0

    Ca = iceberg_constants['form_drag_coefficient_in_air']
    Cw = iceberg_constants['form_drag_coefficient_in_water']
    Cda = iceberg_constants['skin_drag_coefficient_in_air']
    Cdw = iceberg_constants['skin_drag_coefficient_in_water']
    As = iceberg_constants['sail_area']
    Ak = iceberg_constants['keel_area']
    At = iceberg_constants['top_area']
    Ab = iceberg_constants['bottom_area']
    M = iceberg_constants['mass']
    phi = iceberg_constants['latitude']

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
    Fwpx = M * (Amwx + f * Vmwx)
    Fwpy = M * (Amwy - f * Vmwy)

    # Iceberg acceleration
    Ax = (Fax + Fwx + Fcx + Fwpx) / M
    Ay = (Fay + Fwy + Fcy + Fwpy) / M
    
    forces = [Fax, Fay, Fwx, Fwy, Fcx, Fcy, Fwpx, Fwpy] 

    return Ax, Ay, forces
