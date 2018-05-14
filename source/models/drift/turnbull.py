"""Iceberg drift model.

Computes the new x and y components of iceberg position and velocity after one timestep.

Todo:
    * Extend turnbull_drift to include n-layer keel

"""

import numpy as np

def turnbull_drift(iceberg, UA, VA, UW, VW, dt):
    """Computes the new velocity and position of an iceberg after one timestep.
    
    The equations of motion of this function come from the following research paper:
    
    Turnbull, Ian & Fournier, Nicolas & Stolwijk, Michiel & Fosnaes, Tor & Mcgonigal, David. (2014). 
    Operational iceberg drift forecasting in Northwest Greenland. 
    Cold Regions Science and Technology. 
    110. 10.1016/j.coldregions.2014.10.006.
    
    and are best described there.
    
    Notes:
        All units are ANSI. Specifically, velocities are in meters per second, positions are in degrees latitude or 
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
    earth_radius = 6378*1e3  # radius of Earth  (m)
    om = 7.2921e-5  # rotation rate of Earth (rad/s)
    rhoa = 1.225 # density of air (kg/m^3)
    rhow = 1027.5  # density of water (kg/m^3)
    rhoi = 900  # density of iceberg (kg/m^3)
    
    # Iceberg Attributes
    Vx = iceberg.xvels[-1]  # x-component of iceberg velocity (m/s)
    Vy = iceberg.yvels[-1]  # y-component of iceberg velocity (m/s)
    x = iceberg.lons[-1]  # x-component of iceberg position (degrees longitude)
    y = iceberg.lats[-1]  # y-component of iceberg position (degrees latitiude)
    M = iceberg.mass  # iceberg mass (kg)
    Ma = 0.5*M  # added iceberg mass (kg)
    Ca = iceberg.air_drag_coeff  # air drag coefficient
    Cw =  iceberg.water_drag_coeff  # water drag coefficient
    Cdw = iceberg.water_skin_drag_coeff  # skin drag in water coefficient
    Cda = iceberg.air_skin_drag_coeff  # skin drag in air coefficient (Lichey and Hellmer, 2001).
    Ak = iceberg.keel_area  # cross-sectional area of the iceberg's keel (m^2)
    Ab = iceberg.bottom_area  # cross-sectional area of the iceberg's bottom face (m^2)
    As = iceberg.sail_area  # cross-sectional area of the iceberg's keel (m^2) 
    At = iceberg.top_area  # cross-sectional area of the iceberg's top face (m^2)
    
    # Wind and ocean current velocities
    Vax = UA  # u-component of air velocity (m/s)
    Vay = VA  # v-component of air velocity (m/s)
    Vcx = UW  # u-component of ocean current velocity (m/s)
    Vcy = VW  # v-component of ocean current velocity (m/s)
    
    # Air force
    Fax = (0.5*rhoa*Ca*As + rhoa*Cda*At)*abs(Vax - Vx)*(Vax - Vx)  # x-component of the force of wind on iceberg (N)
    Fay = (0.5*rhoa*Ca*As + rhoa*Cda*At)*abs(Vay - Vy)*(Vay - Vy)  # y-component of the force of wind on iceberg (N)
    
    # Water force
    Fwx = (0.5*rhow*Cw*Ak  + rhow*Cdw*Ab)*abs(Vcx - Vx)*(Vcx - Vx)  # x-component of the force of ocean current on iceberg (N)
    Fwy = (0.5*rhow*Cw*Ak + rhow*Cdw*Ab)*abs(Vcy - Vy)*(Vcy - Vy)  # y-component of the force of ocean current on iceberg (N)
      
    # Coriolis force
    f = 2*om*np.sin(np.deg2rad(y))  # Coriolis parameter
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
    ax = (Fax + Fcx + Fwx + Fwpx)/(M + Ma)  # x-component of iceberg acceleration (m/s^2)
    ay = (Fay + Fcy + Fwy + Fwpy)/(M + Ma)  # y-component of iceberg acceleration (m/s^2)
    
    # Iceberg velocity
    Vx_new = Vx + dt*ax  # x-component of iceberg velocity (m/s)
    Vy_new = Vy + dt*ay  # y-component of iceberg velocity (m/s)
    
    # Iceberg position (note the conversion from meters back to degrees)
    y_new = y + dt*Vy_new*(180/(np.pi*earth_radius))  # y-component of iceberg position (degrees latitude)
    x_new = x + dt*Vx_new/(np.cos((((y + y_new)/2)*np.pi)/180))*(180/(np.pi*earth_radius))  # x-component of iceberg position (degrees longitude)
    
    return Vx_new, Vy_new, x_new, y_new
