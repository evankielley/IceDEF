import cmath
import numpy as np
from icedef.constants import *


def newtonian_drift_wrapper(t, lon, lat, vx, vy, **kwargs):

    kwargs['Vcx'] = kwargs.pop('eastward_current').interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values
    kwargs['Vcy'] = kwargs.pop('northward_current').interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values
    kwargs['Vwx'] = kwargs.pop('eastward_wind').interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values
    kwargs['Vwy'] = kwargs.pop('northward_wind').interp(time=t, latitude=lat, longitude=lon, assume_sorted=True).values

    kwargs['phi'] = lat

    ax, ay = newtonian_drift(vx, vy, **kwargs)

    return vx, vy, ax, ay


def newtonian_drift(Vx, Vy, **kwargs):
    """Computes instantaneous iceberg acceleration."""

    # Constants
    Omega = EARTH_ROTATION_RATE
    rhoa = AIR_DENSITY
    rhow = SEAWATER_DENSITY

    Vwx = kwargs.pop('Vwx')
    Vwy = kwargs.pop('Vwy')
    Vcx = kwargs.pop('Vcx')
    Vcy = kwargs.pop('Vcy')

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

    return ax, ay


def dev_newtonian_drift(iceberg_velocity, current_velocity, wind_velocity, **kwargs):
    """Computes instantaneous iceberg acceleration."""

    # Override forces with kwargs
    zero_wind_force = kwargs.pop('zero_wind_force', False)
    zero_current_force = kwargs.pop('zero_current_force', False)
    zero_coriolis_force = kwargs.pop('zero_coriolis_force', False)
    zero_pressure_force = kwargs.pop('zero_pressure_force', False)

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
    if zero_wind_force:
        Fax, Fay = 0, 0
    else:
        Fax = (0.5 * rhoa * Ca * As + rhoa * Cda * At) * np.sqrt((Vwx - Vx)**2 + (Vwy - Vy)**2) * (Vwx - Vx)
        Fay = (0.5 * rhoa * Ca * As + rhoa * Cda * At) * np.sqrt((Vwx - Vx)**2 + (Vwy - Vy)**2) * (Vwy - Vy)

    # Current force

    ekman = kwargs.pop('ekman', False)

    if zero_current_force:
        Fwx, Fwy = 0, 0

    elif ekman:

        Fwx_list = []
        Fwy_list = []
        Vcx_list = []
        Vcy_list = []

        depth_vec = kwargs.pop('depth_vec', np.arange(0, -110, -10))

        u_vec, v_vec = compute_ekman_spiral(wind_velocity, current_velocity, depth_vec)

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
    if zero_coriolis_force:
        Fcx, Fcy = 0, 0
    else:
        Fcx = f * M * Vy
        Fcy = -f * M * Vx

    # Water pressure force
    if zero_pressure_force:
        Fwpx, Fwpy = 0, 0
    else:
        Vmwx = Vcx
        Vmwy = Vcy
        Amwx = Amwx
        Amwy = Amwy
        Fwpx = M * (Amwx - f * Vmwy)
        Fwpy = M * (Amwy + f * Vmwx)

    # Iceberg acceleration
    Ax = (Fax + Fwx + Fcx + Fwpx) / (M + 0.5 * M)
    Ay = (Fay + Fwy + Fcy + Fwpy) / (M + 0.5 * M)

    forces = [Fax, Fay, Fwx, Fwy, Fcx, Fcy, Fwpx, Fwpy]

    return Ax, Ay, forces


def compute_ekman_velocity(wind, depth, latitude=50):
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


def analytical_drift(iceberg_position, current_velocity, wind_velocity, **kwargs):
    y, x = iceberg_position
    vau, vav = wind_velocity
    vwu, vwv = current_velocity
    l = kwargs.pop('iceberg_length', 160)
    w = kwargs.pop('iceberg_width', 160)
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
