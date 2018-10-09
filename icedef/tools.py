#import numpy as np
#from icedef.constants import *
#from scipy.interpolate import interp1d


def dx_to_dlon(dx, lat):
    return dx / (np.cos((lat * np.pi) / 180)) * (180 / (np.pi * EARTH_RADIUS))


def dy_to_dlat(dy):
    return dy * (180 / (np.pi * EARTH_RADIUS))


def dlon_to_dx(dlon, lat):
    return dlon * np.cos((lat * np.pi) / 180) * (np.pi * EARTH_RADIUS) / 180


def dlat_to_dy(dlat):
    return dlat * (np.pi * EARTH_RADIUS) / 180


def deg2m(lon0, lonn, lat0, latn):
    """ Calculates x and y distances in meters.
    """

    # constants
    PI = np.pi
    Ce = 40075160  # circumference of Earth around equator
    Cp = 40008000  # circumference of Earth around poles
    Clat0 = Ce * np.cos(lat0*PI/180)  # circumference of Earth at lat0

    dlat = latn - lat0
    dlon = lonn - lon0
    dx = dlon * Clat0 / 360
    dy = dlat * Cp / 360

    return dx, dy


def compute_mse(simulation_point, observation_vectors, reference_time):

    sim_x, sim_y, sim_t = simulation_point
    obs_x_vec, obs_y_vec, obs_t_vec = observation_vectors

    t0 = reference_time

    sim_t = (sim_t - t0).item().total_seconds()
    obs_t_vec = [(obs_t - t0).item().total_seconds() for obs_t in obs_t_vec]

    obs_x_interpolant = interp1d(obs_t_vec, obs_x_vec)
    obs_y_interpolant = interp1d(obs_t_vec, obs_y_vec)

    obs_x = obs_x_interpolant(sim_t)
    obs_y = obs_y_interpolant(sim_t)

    mse = np.sqrt((obs_x - sim_x)**2 + (obs_y - sim_y)**2)

    return mse
