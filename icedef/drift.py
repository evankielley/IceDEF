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

    ekman = kwargs.pop('ekman', False)

    if ekman:

        Fwx_list = []
        Fwy_list = []

        depth_vec = kwargs.pop('depth_vec', np.arange(0, -110, -10))

        u_vec, v_vec = compute_ekman_spiral(wind_velocity, current_velocity, depth_vec)

        for i in range(len(u_vec)):

            Vcx, Vcy = u_vec[i], v_vec[i]

            Fwx_list.append((0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * np.sqrt((Vcx - Vx)**2 + (Vcy - Vy)**2) * (Vcx - Vx))
            Fwy_list.append((0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * np.sqrt((Vcx - Vx)**2 + (Vcy - Vy)**2) * (Vcy - Vy))

        Fwx = np.mean(np.array(Fwx_list))
        Fwy = np.mean(np.array(Fwy_list))

    else:

        Fwx = (0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * np.sqrt((Vcx - Vx) ** 2 + (Vcy - Vy) ** 2) * (Vcx - Vx)
        Fwy = (0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * np.sqrt((Vcx - Vx) ** 2 + (Vcy - Vy) ** 2) * (Vcy - Vy)

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


def compute_ekman_velocity(wind, depth, latitude=50):

    u_wind, v_wind = wind
    z = depth
    phi = latitude

    rho_air = 1.225  # density of air (kg/m^3)
    rho_water = 1028  # density of seawater (kg/m^3)
    Cd =  1.3e-3 # ranges from (1.1 - 1.5) x 10^-3
    Omega = 7.2910e-5  # rotation rate of Earth (s^-1)
    Az = 5e-2   # m^2 s^−1;

    f = lambda phi : 2 * Omega * np.sin(phi)  # Coriolis parameter

    tau_x = lambda U, V : rho_air * Cd * U * np.sqrt(U**2 + V**2)  # wind stress x-component
    tau_y = lambda U, V : rho_air * Cd * V * np.sqrt(U**2 + V**2)  # wind stress y-component

    V0x = tau_x(u_wind, v_wind) / np.sqrt(rho_water**2 * np.abs(f(phi)) * Az)
    V0y = tau_y(u_wind, v_wind) / np.sqrt(rho_water**2 * np.abs(f(phi)) * Az)
    V0 = np.sqrt(V0x**2 + V0y**2)
    theta = np.pi/2 - np.arctan2(V0y, V0x)

    a = np.sqrt(abs(f(phi))/(2 * Az))

    # note: clockwise rotation
    u_ekman = V0 * np.exp(a * z) * np.cos(np.pi / 4 + a * z) * (np.cos(theta) + np.sin(theta))
    v_ekman = V0 * np.exp(a * z) * np.sin(np.pi / 4 + a * z) * (np.cos(theta) - np.sin(theta))

    return u_ekman, v_ekman


def compute_ekman_spiral(wind, surface_current, depth_vec):

    u_ekman_vec = np.zeros(len(depth_vec))
    v_ekman_vec = np.zeros(len(depth_vec))

    for i, depth in enumerate(depth_vec):

        u_ekman_vec[i], v_ekman_vec[i] = compute_ekman_velocity(wind, depth)

    u_surface_current, v_surface_current = surface_current

    u_barotropic = u_surface_current - u_ekman_vec[0]
    v_barotropic = v_surface_current - v_ekman_vec[0]

    u_current_vec = u_barotropic + u_ekman_vec
    v_current_vec = v_barotropic + v_ekman_vec

    return u_current_vec, v_current_vec
