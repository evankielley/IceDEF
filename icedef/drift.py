"""Iceberg drift functions"""

from icedef.constants import *
import numpy as np


def newtonian_drift(velocity, **kwargs):
    """Computes instantaneous iceberg acceleration."""

    # Constants
    om = EARTH_ROTATION_RATE
    rhoa = AIR_DENSITY
    rhow = SEAWATER_DENSITY
    Ca = ICEBERG_FORM_DRAG_COEFFICIENT_IN_AIR
    Cw = ICEBERG_FORM_DRAG_COEFFICIENT_IN_WATER
    Cda = ICEBERG_SKIN_DRAG_COEFFICIENT_IN_AIR
    Cdw = ICEBERG_SKIN_DRAG_COEFFICIENT_IN_WATER

    # Kwargs
    Vwx, Vwy = kwargs.get('wind_velocity')
    Vcx, Vcy = kwargs.get('current_velocity')
    As = kwargs.get('iceberg_sail_area')
    Ak = kwargs.get('iceberg_keel_area')
    At = kwargs.get('iceberg_top_area')
    Ab = kwargs.get('iceberg_bottom_area')
    M = kwargs.get('iceberg_mass')
    phi = kwargs.get('iceberg_latitude')

    # Args
    Vx = velocity[0]
    Vy = velocity[1]

    # Wind force
    Fax = 0.5 * rhoa * Ca * As + rhoa * Cda * At * abs(Vwx - Vx) * (Vwx - Vx)
    Fay = 0.5 * rhoa * Ca * As + rhoa * Cda * At * abs(Vwy - Vy) * (Vwy - Vy)

    # Current force
    Fwx = 0.5 * rhow * Cw * Ak + rhow * Cdw * Ab * abs(Vcx - Vx) * (Vcx - Vx)
    Fwy = 0.5 * rhow * Cw * Ak + rhow * Cdw * Ab * abs(Vcy - Vy) * (Vcy - Vy)

    # Coriolis force
    f = 2 * om * np.sin(np.deg2rad(phi))
    Fcx = f * M * Vy
    Fcy = -f * M * Vx

    # Water pressure force
    Vmwx = Vcx
    Vmwy = Vcy
    Amwx = 0
    Amwy = 0
    Fwpx = M * (Amwx + f * Vmwx)
    Fwpy = M * (Amwy + f * Vmwy)

    # Iceberg acceleration
    Ax = (Fax + Fwx + Fcx + Fwpx) / M
    Ay = (Fay + Fwy + Fcy + Fwpy) / M

    return Ax, Ay
