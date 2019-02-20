"""This module contains iceberg drift models.
"""

import cmath
import numpy as np
from icedef.constants import *


def newtonian_drift_wrapper(t, lon, lat, vx, vy, **kwargs):
    """This function performs interpolations for current and wind velocities and then runs the drift model.

    Args:
        t (numpy.datetime64): time.
        lon (float): longitude.
        lat (float): latitude.
        vx (float): x-component of iceberg velocity in m/s.
        vy (float): y-component of iceberg velocity in m/s.
        **kwargs: coming soon - see source code for now.

    Returns:
        vx (float): new x-component of iceberg velocity in m/s.
        vy (float): new y-component of iceberg velocity in m/s.
        ax (float): x-component of iceberg acceleration in m/s.
        ay (float): y-component of iceberg acceleration in m/s.
    """

    dt = kwargs.pop('time_step', np.timedelta64(300, 's'))

    fast_interpolation = kwargs.pop('fast_interpolation', True)

    if fast_interpolation:

        current_interpolator = kwargs.pop('current_interpolator')
        Vcx, Vcy = current_interpolator((t, lat, lon))
        wind_interpolator = kwargs.pop('wind_interpolator')
        Vwx, Vwy = wind_interpolator((t, lat, lon))

        Vcx_left, Vcy_left = current_interpolator((t - dt, lat, lon))
        Vcx_right, Vcy_right = current_interpolator((t + dt, lat, lon))

    else:

        Vcxs = kwargs.pop('eastward_current')
        Vcys = kwargs.pop('northward_current')
        Vwxs = kwargs.pop('eastward_wind')
        Vwys = kwargs.pop('northward_wind')

        Vcx = Vcxs.interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values
        Vcy = Vcys.interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values
        Vwx = Vwxs.interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values
        Vwy = Vwys.interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values

        Vcx_left = Vcxs.interp(time=t - dt, latitude=lat, longitude=lon, assume_sorted=True).values
        Vcx_right = Vcxs.interp(time=t + dt, latitude=lat, longitude=lon, assume_sorted=True).values
        Vcy_left = Vcys.interp(time=t - dt, latitude=lat, longitude=lon, assume_sorted=True).values
        Vcy_right = Vcys.interp(time=t + dt, latitude=lat, longitude=lon, assume_sorted=True).values

    Amwx = (Vcx_right - Vcx_left) / (dt.item().total_seconds() * 2)
    Amwy = (Vcy_right - Vcy_left) / (dt.item().total_seconds() * 2)

    kwargs['Vcx'] = Vcx
    kwargs['Vcy'] = Vcy
    kwargs['Vwx'] = Vwx
    kwargs['Vwy'] = Vwy
    kwargs['Amwx'] = Amwx
    kwargs['Amwy'] = Amwy

    kwargs['phi'] = lat

    ax, ay = newtonian_drift(vx, vy, **kwargs)

    return vx, vy, ax, ay


def newtonian_drift(Vx, Vy, **kwargs):
    """This function computes iceberg acceleration using a general Newtonian drift model.

    Args:
        Vx (float): x-component of iceberg velocity in m/s.
        Vy (float): y-component of iceberg velocity in m/s.
        **kwargs: coming soon - see source code for now.

    Returns:
        ax (float): x-component of iceberg acceleration in m/s.
        ay (float): y-component of iceberg acceleration in m/s.
    """

    # Constants
    Omega = EARTH_ROTATION_RATE
    rhoa = AIR_DENSITY
    rhow = SEAWATER_DENSITY

    Vwx = kwargs.pop('Vwx')
    Vwy = kwargs.pop('Vwy')
    Vcx = kwargs.pop('Vcx')
    Vcy = kwargs.pop('Vcy')

    Amwx = kwargs.pop('Amwx', 0)
    Amwy = kwargs.pop('Amwy', 0)

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
    ekman = kwargs.pop('ekman', False)

    if ekman:

        Fwx_list = []
        Fwy_list = []
        Vcx_list = []
        Vcy_list = []

        depth_vec = kwargs.pop('depth_vec', np.arange(0, -110, -10))

        u_vec, v_vec = compute_ekman_spiral((Vwx, Vwy), (Vcx, Vcy), depth_vec)

        for i in range(len(u_vec)):

            Vcx, Vcy = u_vec[i], v_vec[i]

            Vcx_list.append(Vcx)
            Vcy_list.append(Vcy)

            Fwx = (0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * np.sqrt((Vcx - Vx)**2 + (Vcy - Vy)**2) * (Vcx - Vx)
            Fwy = (0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * np.sqrt((Vcx - Vx)**2 + (Vcy - Vy)**2) * (Vcy - Vy)

            Fwx_list.append(Fwx)
            Fwy_list.append(Fwy)

        Vcx = np.mean(np.array(Vcx_list))
        Vcy = np.mean(np.array(Vcy_list))

        Fwx = np.mean(np.array(Fwx_list))
        Fwy = np.mean(np.array(Fwy_list))

    else:
        Fwx = (0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * np.sqrt((Vcx - Vx)**2 + (Vcy - Vy)**2) * (Vcx - Vx)
        Fwy = (0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * np.sqrt((Vcx - Vx)**2 + (Vcy - Vy)**2) * (Vcy - Vy)

    # Coriolis Parameter
    f = 2 * Omega * np.sin(np.deg2rad(phi))

    # Coriolis force
    Fcx = f * M * Vy
    Fcy = -f * M * Vx

    # Water Pressure Gradient Force
    Vmwx = Vcx
    Vmwy = Vcy
    Amwx = Amwx
    Amwy = Amwy
    Fwpx = M * (Amwx - f * Vmwy)
    Fwpy = M * (Amwy + f * Vmwx)

    # Iceberg acceleration
    ax = (Fax + Fwx + Fcx + Fwpx) / (M + 0.5 * M)
    ay = (Fay + Fwy + Fcy + Fwpy) / (M + 0.5 * M)

    log = kwargs.pop('log', None)

    if log is not None:
        log.info('{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}'.format(
            Fax, Fay, Fwx, Fwy, Fcx, Fcy, Fwpx, Fwpy))

    return ax, ay


def compute_ekman_velocity(wind, depth, latitude=50):
    """This function computes Ekman velocity at some depth.

    Args:
        wind (tuple of float): components (x, y) of wind velocity in m/s.
        depth (float): depth below the sea surface (down is negative) in m.
        latitude: latitude.

    Returns:
        u_ekman (float): x-component of Ekman velocity at specified depth.
        v_ekman (float): y-component of Ekman velocity at specified depth.
    """

    u_wind, v_wind = wind
    z = depth
    phi = latitude

    rho_air = 1.225  # density of air (kg/m^3)
    rho_water = 1028  # density of seawater (kg/m^3)
    Cd = 1.3e-3  # ranges from (1.1 - 1.5) x 10^-3
    Omega = 7.2910e-5  # rotation rate of Earth (s^-1)
    Az = 5e-2  # m^2 s^−1;

    f = lambda phi: 2 * Omega * np.sin(phi)  # Coriolis parameter

    tau_x = lambda U, V: rho_air * Cd * U * np.sqrt(U ** 2 + V ** 2)  # wind stress x-component
    tau_y = lambda U, V: rho_air * Cd * V * np.sqrt(U ** 2 + V ** 2)  # wind stress y-component

    V0x = tau_x(u_wind, v_wind) / np.sqrt(rho_water ** 2 * np.abs(f(phi)) * Az)
    V0y = tau_y(u_wind, v_wind) / np.sqrt(rho_water ** 2 * np.abs(f(phi)) * Az)
    V0 = np.sqrt(V0x ** 2 + V0y ** 2)
    theta = np.pi / 2 - np.arctan2(V0y, V0x)

    a = np.sqrt(abs(f(phi)) / (2 * Az))

    # note: clockwise rotation
    u_ekman = V0 * np.exp(a * z) * np.cos(np.pi / 4 + a * z) * (np.cos(theta) + np.sin(theta))
    v_ekman = V0 * np.exp(a * z) * np.sin(np.pi / 4 + a * z) * (np.cos(theta) - np.sin(theta))

    return u_ekman, v_ekman


def compute_ekman_spiral(wind, surface_current, depth_vec, latitude=50):
    """This function computes an Ekman spiral.

    Args:
        wind (tuple of float): components (x, y) of wind velocity in m/s.
        surface_current (tuple of float): components (x, y) of current velocity at the surface in m/s.
        depth_vec (list of float): depths in m to compute Ekman velocity at.
        latitude: latitude.

    Returns:
        u_current_vec (list of float): x-components of Ekman velocity at all depths specified.
        v_current_vec (list of float): y-components of Ekman velocity at all depths specified.
    """

    u_ekman_vec = np.zeros(len(depth_vec))
    v_ekman_vec = np.zeros(len(depth_vec))

    u_ekman_surface, v_ekman_surface = compute_ekman_velocity(wind, 0, latitude)

    u_surface_current, v_surface_current = surface_current

    u_barotropic = u_surface_current - u_ekman_surface
    v_barotropic = v_surface_current - v_ekman_surface

    for i, depth in enumerate(depth_vec):
        u_ekman_vec[i], v_ekman_vec[i] = compute_ekman_velocity(wind, depth, latitude)

    u_current_vec = u_barotropic + u_ekman_vec
    v_current_vec = v_barotropic + v_ekman_vec

    return u_current_vec, v_current_vec


def analytical_drift_wrapper(t, lon, lat, **kwargs):
    """This function performs interpolations for current and wind velocities and then runs the drift model.

    Args:
        t (numpy.datetime64): time.
        lon (float): longitude.
        lat (float): latitude.
        vx (float): x-component of iceberg velocity in m/s.
        vy (float): y-component of iceberg velocity in m/s.
        **kwargs: coming soon - see source code for now.

    Returns:
        vx (float): new x-component of iceberg velocity in m/s.
        vy (float): new y-component of iceberg velocity in m/s.

    """

    fast_interpolation = kwargs.pop('fast_interpolation', True)

    if fast_interpolation:

        current_interpolator = kwargs.pop('current_interpolator')
        wind_interpolator = kwargs.pop('wind_interpolator')

        vwu, vwv = current_interpolator((t, lat, lon))
        vau, vav = wind_interpolator((t, lat, lon))

    else:

        vwus = kwargs.pop('eastward_current')
        vwvs = kwargs.pop('northward_current')
        vaus = kwargs.pop('eastward_wind')
        vavs = kwargs.pop('northward_wind')

        vwu = vwus.interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values
        vwv = vwvs.interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values
        vau = vaus.interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values
        vav = vavs.interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values

    kwargs['vwu'] = vwu
    kwargs['vwv'] = vwv
    kwargs['vau'] = vau
    kwargs['vav'] = vav

    vx, vy = analytical_drift(lon, lat, **kwargs)

    return vx, vy


def analytical_drift(x, y, **kwargs):
    """This function computes the velocity of an iceberg using an analytical drift model.

    Args:
        x (float): longitude
        y (float): latitude
        **kwargs: coming soon - see source code for now.

    Returns:
        viu (float): x-component of iceberg velocity in m/s.
        viv (float): y-component of iceberg velocity in m/s.
    """

    vwu = kwargs.pop('vwu')
    vwv = kwargs.pop('vwv')
    vau = kwargs.pop('vau')
    vav = kwargs.pop('vav')

    l = kwargs.pop('waterline_length', 160)
    w = kwargs.pop('waterline_length', 160)  # note: currently all icebergs are cuboid

    Cw = kwargs.pop('form_drag_coefficient_in_water', 0.9)
    Ca = kwargs.pop('form_drag_coefficient_in_air', 1.3)

    Omega = EARTH_ROTATION_RATE
    rhow = SEAWATER_DENSITY
    rhoa = AIR_DENSITY
    rhoi = ICEBERG_DENSITY

    gamma = np.sqrt(rhoa * (rhow - rhoi) / rhow / rhoi * (Ca / Cw))
    S = np.pi * ((l * w) / (l + w))
    f = 2 * Omega * np.sin((np.abs(y) * np.pi) / 180)
    Lambda = np.sqrt(2) * Cw * (gamma * np.sqrt(vau ** 2 + vav ** 2)) / (f * S)

    if Lambda < 0.1:
        alpha = Lambda * (Lambda**4 * (Lambda**4 * (Lambda**4 * (-0.0386699020961393 * Lambda**4 +
                0.055242717280199) - 0.0883883476483184) + 0.176776695296637) - 0.707106781186548)

    else:
        alpha = np.multiply(np.divide(np.sqrt(2), np.power(Lambda, 3)), (1 - np.sqrt(1 + np.power(Lambda, 4))))

    if Lambda < 0.6:
        beta = Lambda**3 * (Lambda**4 * (Lambda**4 * (Lambda**4 * (Lambda**4 *
                (Lambda**4 * (Lambda**4 * (Lambda**4 * (Lambda**4 * (0.0153268598203613 *
                Lambda**4 - 0.0151656272365985) + 0.0180267866272764) + 0.0219176256311202) -
                0.0274446790511418) + 0.0357675015202851) - 0.0493731785691779) + 0.0745776683282687) -
                0.132582521472478) + 0.353553390593274)

    else:
        beta = np.real(np.multiply(np.divide(1, np.power(Lambda, 3)), cmath.sqrt(np.multiply((4 +
                np.power(Lambda, 4)), cmath.sqrt(1 + np.power(Lambda, 4))) - 3 * np.power(Lambda, 4) - 4)))

    viu = vwu + gamma * (-alpha * vav + beta * vau)
    viv = vwv + gamma * (alpha * vau + beta * vav)

    return viu, viv
