import numpy as np

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