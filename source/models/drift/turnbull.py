import numpy as np

def turnbull_drift(iceberg, UA, VA, UW, VW, dt):
    
    earth_radius = 6378*1e3
    om = 7.2921e-5  # rotation rate of earth in rad/s
    rhoa = 1.225 # density of air (kg/m^3)
    rhow = 1027.5  # density of water (kg/m^3)
    rhoi = 900


    Vax = UA
    Vay = VA
    Vcx = UW
    Vcy = VW
        
    Vx = iceberg.xvels[-1]
    Vy = iceberg.yvels[-1]
    x = iceberg.lons[-1]
    y = iceberg.lats[-1]
    l = iceberg.length[-1]
    w = iceberg.width[-1]
    h = iceberg.height[-1]
    berg_mass = iceberg.mass
    Ma = 0.5*berg_mass  # added mass
        

    # calculate air density
    # TODO: effect of temperature and humidity

    rhoi = iceberg.density
    Ca = iceberg.air_drag_coeff  # air drag coefficient
    Cw =  iceberg.water_drag_coeff  # water drag coefficient
    
    # Keel info
    Cdw = iceberg.water_skin_drag_coeff  # skin drag in water coefficient
    keel_shape = iceberg.keel_shape
    #Ak = iceberg.keel_area
    #Ab = iceberg.bottom_area
    Ak = (rhoi/rhow)*(2/np.pi)*(l+w)*h / 2
    Ab = 0
    
    # Sail info
    Cda = iceberg.air_skin_drag_coeff  # skin drag in air coefficient (Lichey and Hellmer, 2001).
    sail_shape = iceberg.sail_shape
    #As = iceberg.sail_area  
    #At = iceberg.top_area
    As = ((rhow - rhoi)/rhoi)*Ak / 2
    At = l*w
    
    
    # Air force
    Fax = (0.5*rhoa*Ca*As + rhoa*Cda*At)*abs(Vax - Vx)*(Vax - Vx)
    Fay = (0.5*rhoa*Ca*As + rhoa*Cda*At)*abs(Vay - Vy)*(Vay - Vy)
    

    # Water force
    # TODO: extend to n-layer keel model
    Fwx = (0.5*rhow*Cw*Ak  + rhow*Cdw*Ab)*abs(Vcx - Vx)*(Vcx - Vx) 
    Fwy = (0.5*rhow*Cw*Ak + rhow*Cdw*Ab)*abs(Vcy - Vy)*(Vcy - Vy)
    
    
    # Coriolis force
    f = 2*om*np.sin(np.deg2rad(y))  # Coriolis parameter (y is latitude in degrees)
    Fcx = f*Vy*berg_mass
    Fcy = -f*Vx*berg_mass
    
    
    # Water pressure gradient
    # Mean water current down to the iceberg keel
    Vwmx = 0
    Vwmy = 0
    
    # acceleration (time-derivative) of Vmw
    Amwx = 0
    Amwy = 0
        
    Fwpx = berg_mass*(Amwx + f*Vwmx)
    Fwpy = berg_mass*(Amwy - f*Vwmy)
    
    
    # Iceberg acceleration                
    ax = (Fax + Fcx + Fwx + Fwpx)/(berg_mass + Ma)
    ay = (Fay + Fcy + Fwy + Fwpy)/(berg_mass + Ma)
    
    # Iceberg velocity
    Vx += dt*ax
    Vy += dt*ay
    
    # Iceberg position
    y_new = y + dt*Vy*(180/(np.pi*earth_radius))
    x_new = x + dt*Vx/(np.cos((((y + y_new)/2)*np.pi)/180))*(180/(np.pi*earth_radius))
    
    return Vx, Vy, x_new, y_new
