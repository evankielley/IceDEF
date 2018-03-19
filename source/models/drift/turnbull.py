import numpy as np

def turnbull_drift(t, dt, vx, vy, x, y, l, w, h, UA, VA, UW, VW):

    Vax = UA
    Vay = VA
    Vcx = UW
    Vcy = VW
        
    Vx = vx
    Vy = vy
        

    # calculate air density
    # TODO: effect of temperature and humidity
    rhoa = 1.225 # density of air (kg/m^3)
    rhow = 1027.5  # density of water (kg/m^3)
    rhoi = 900
    Ca = 1.2e-3 #+ 0.1)/2  # air drag coefficient
    Cw =  0.5 #(5.5e-3+ 2.5)/2    # water drag coefficient
    
    # Keel info
    Cdw = 5.0e-4  # skin drag in water coefficient
    keel_shape = 'triangular'
    #keel_shape = 'rectangular'
    Ab = 0

    if keel_shape is 'triangular':
        Ak = (rhoi/rhow)*(2/np.pi)*(l+w)*h / 2  # from Wagner   
    elif keel_shape is 'rectangular':
        Ak = (rhoi/rhow)*(2/np.pi)*(l+w)*h  # from Wagner      
    else:
        print('Invalid keel shape')
        
    
    # Sail info
    Cda = 2.5e-4  # skin drag in air coefficient (Lichey and Hellmer, 2001).
    sail_shape = 'rectangular'
    As = ((rhow - rhoi)/rhoi)*Ak / 2   # from Wagner   
    At = l*w  # area of top of sail
    
    om = 7.2921e-5  # rotation rate of earth in rad/s
    f = 2*om*np.sin(np.deg2rad(y))  # Coriolis parameter (y is latitude in degrees)

    berg_mass = l*w*h*rhoi
    Ma = 0.5*berg_mass  # added mass
    
    earth_radius = 6378*1e3

    
    # Air force
    Fax = (0.5*rhoa*Ca*As + rhoa*Cda*At)*abs(Vax - Vx)*(Vax - Vx)
    Fay = (0.5*rhoa*Ca*As + rhoa*Cda*At)*abs(Vay - Vy)*(Vay - Vy)
    

    # Water force
    # TODO: extend to n-layer keel model
    Fwx = (0.5*rhow*Cw*Ak  + rhow*Cdw*Ab)*abs(Vcx - Vx)*(Vcx - Vx) 
    Fwy = (0.5*rhow*Cw*Ak + rhow*Cdw*Ab)*abs(Vcy - Vy)*(Vcy - Vy)
    
    
    # Coriolis force
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
    vx += dt*ax
    vy += dt*ay
    
    # Iceberg position
    y_new = y + dt*vy*(180/(np.pi*earth_radius))
    x_new = x + dt*vx/(np.cos((((y + y_new)/2)*np.pi)/180))*(180/(np.pi*earth_radius))
    
    return vx, vy, x_new, y_new
