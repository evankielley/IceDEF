"""Iceberg drift model.

Computes the new x and y components of iceberg position and velocity after one timestep.

Todo:
    * Extend turnbull_drift to include n-layer keel

"""

import numpy as np
import netCDF4 as nc

def drift(iceberg, Vax, Vay, Vcx, Vcy, dt):
    """Computes the new velocity and position of an iceberg after one timestep.
    
    The equations of motion of this function come from the following research paper:
    
    Turnbull, Ian & Fournier, Nicolas & Stolwijk, Michiel & Fosnaes, Tor & Mcgonigal, David. (2014). 
    Operational iceberg drift forecasting in Northwest Greenland. 
    Cold Regions Science and Technology. 
    110. 10.1016/j.coldregions.2014.10.006.
    
    and are best described there.
    
    Notes:
        All units are SI. Specifically, velocities are in meters per second, positions are in degrees latitude or 
        longitude, size dimensions are in meters, masses are in kilograms, times are in seconds, and forces are in Newtons.
    
    Args:
        iceberg (obj): Contains variables that describe the iceberg's velocity, position, size dimensions, and mass.
        UA (float): U-component of air velocity.
        VA (float): V-component of air velocity.
        UW (float): U-component of water velocity.
        VW (float): V-component of water velocity.
        dt (float): Timestep size.
        
    Returns:
        Vx_new (float): New x-component of iceberg velocity
        Vy_new (float): New y-component of iceberg velocity
        x_new (float): New x-component of iceberg position.
        y_new (float): New y-component of iceberg position.
    
    """
    
    # Physical constants
    om = 7.2921e-5  # rotation rate of Earth (rad/s)
    rhoa = 1.225 # density of air (kg/m^3)
    rhow = 1027.5  # density of water (kg/m^3)
    rhoi = iceberg.rho  # density of iceberg (kg/m^3)
    
    # Iceberg Attributes
    T = iceberg.T  # time of the iceberg (datetime)
    X = iceberg.X  # x-component of iceberg position (degrees longitude)
    Y = iceberg.Y  # y-component of iceberg position (degrees latitiude)
    Vx = iceberg.Vx  # x-component of iceberg velocity (m/s)
    Vy = iceberg.Vy  # y-component of iceberg velocity (m/s)
    M = iceberg.M  # iceberg mass (kg)
    Ma = 0.5*M  # added iceberg mass (kg)
    Cda = iceberg.Cda  # air drag coefficient
    Cdw =  iceberg.Cdw  # water drag coefficient
    Csdw = iceberg.Csdw  # skin drag in water coefficient
    Csda = iceberg.Csda  # skin drag in air coefficient (Lichey and Hellmer, 2001).
    Ak = iceberg.Ak  # cross-sectional area of the iceberg's keel (m^2)
    Ab = iceberg.Ab  # cross-sectional area of the iceberg's bottom face (m^2)
    As = iceberg.As  # cross-sectional area of the iceberg's keel (m^2) 
    At = iceberg.At  # cross-sectional area of the iceberg's top face (m^2)


    # Air force
    Fax = (0.5*rhoa*Cda*As + rhoa*Csda*At)*abs(Vax - Vx)*(Vax - Vx)  # x-component of the force of wind on iceberg (N)
    Fay = (0.5*rhoa*Cda*As + rhoa*Csda*At)*abs(Vay - Vy)*(Vay - Vy)  # y-component of the force of wind on iceberg (N)
    
    # Water force
    Fwx = (0.5*rhow*Cdw*Ak  + rhow*Csdw*Ab)*abs(Vcx - Vx)*(Vcx - Vx)  # x-component of the force of ocean current on iceberg (N)
    Fwy = (0.5*rhow*Cdw*Ak + rhow*Csdw*Ab)*abs(Vcy - Vy)*(Vcy - Vy)  # y-component of the force of ocean current on iceberg (N)
      
    # Coriolis force
    f = 2*om*np.sin(np.deg2rad(Y))  # Coriolis parameter
    Fcx = f*Vy*M  # x-component of the Coriolis force on iceberg
    Fcy = -f*Vx*M  # y-component of the Coriolis force on iceberg
      
    # Water pressure gradient force
    Vwmx = 0  # x-component of mean water current down to the iceberg keel (m/s)
    Vwmy = 0  # y-component of mean water current down to the iceberg keel (m/s)
    Amwx = 0  # x-component of acceleration (time-derivative) of Vmw (m/s^2)
    Amwy = 0  # y-component of acceleration (time-derivative) of Vmw (m/s^2)    
    Fwpx = M*(Amwx + f*Vwmx)  # x-component of water pressure gradient force (N)
    Fwpy = M*(Amwy - f*Vwmy)  # y-component of water pressure gradient force (N)
      
    # Iceberg acceleration                
    Ax = (Fax + Fcx + Fwx + Fwpx)/(M + Ma)  # x-component of iceberg acceleration (m/s^2)
    Ay = (Fay + Fcy + Fwy + Fwpy)/(M + Ma)  # y-component of iceberg acceleration (m/s^2)
    
    # Iceberg velocity
    Vx_new = Vx + dt*Ax  # x-component of iceberg velocity (m/s)
    Vy_new = Vy + dt*Ay  # y-component of iceberg velocity (m/s)
    

    return Vx_new, Vy_new
    
    
