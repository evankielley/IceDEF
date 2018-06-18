from datetime import timedelta
import netCDF4 as nc
import numpy as np


def get_interpolated_values(berg, ocean_data, atm_data):
    # convert time to format used in metocean fields
    T_ocean = nc.date2num(berg.T, ocean_data.t_units, ocean_data.t_calendar)
    T_atm = nc.date2num(berg.T, atm_data.t_units, atm_data.t_calendar)

    # interpolate current and wind velocities for berg's [t, y, x]
    Vcx = ocean_data.iUW([T_ocean, berg.Y, berg.X])[0] 
    Vcy = ocean_data.iVW([T_ocean, berg.Y, berg.X])[0] 
    Vax = atm_data.iUA([T_atm, berg.Y, berg.X])[0]  
    Vay = atm_data.iVA([T_atm, berg.Y, berg.X])[0]
    
    if abs(Vcx) > 100 or abs(Vcy) > 100:
        print(f'Beware, current speeds are ({Vcx},{Vcy})')
        
    elif abs(Vax) > 200 or abs(Vay) > 200:
        print(f'Beware, wind speeds are ({Vax},{Vay})')
        
    return Vcx, Vcy, Vax, Vay



def get_bounding_box(ocean_data, atm_data):
    drift_xmin = max(min(ocean_data.lons), min(atm_data.lons))
    drift_xmax = min(max(ocean_data.lons), max(atm_data.lons))
    drift_ymin = max(min(ocean_data.lats), min(atm_data.lats))
    drift_ymax = min(max(ocean_data.lats), max(atm_data.lats))
    x_bounds = [drift_xmin, drift_xmax]
    y_bounds = [drift_ymin, drift_ymax]
    return x_bounds, y_bounds



def get_all_timesteps(t0, tf, dt):
    """This function gets a list of all datetimes between t0 and tf in increments of dt
    
    Args: 
        t0 (datetime.datetime): start time
        tf (datetime.datetime): end time
        dt (float): time step (hours)
        
    Returns:
        t_all (list of datetime.datetime): all times between t0 and tf, incremented by dt 
    """
    
    tspan = tf - t0  # timedelta
    tspan_hr = tspan.days*24 + tspan.seconds/3600
    dt_s = timedelta(hours = dt)  # dt in seconds as timedelta
    num_steps = int(-(-(tspan_hr+dt)//dt)) # negatives give ceiling
    
    t_all = [t0 + i*dt_s for i in range(num_steps)]
    
    return t_all



def euler(berg, ocean_data, atm_data, drift, t_all):
    
    Re = 6378*1e3  # radius of Earth
    
    # determine bounding box for drift simulation
    x_bounds, y_bounds = get_bounding_box(ocean_data, atm_data)
    
    # get timestep
    tdelta = t_all[1] - t_all[0]
    dt = tdelta.seconds
    
    # drift iceberg for all t in t_all
    for t in t_all:
        
        berg.T = t
    
        Vcx, Vcy, Vax, Vay = get_interpolated_values(berg, ocean_data, atm_data)
        
        
        # integrate
        # explicit Euler forward scheme
        berg.Vx += dt*berg.Ax
        berg.Vy += dt*berg.Ay   
        berg.Ax, berg.Ay = drift(berg, Vax, Vay, Vcx, Vcy)

         
        # calculate position in degrees lat/lon
        y_ = berg.Y
        berg.Y += dt*berg.Vy*(180/(np.pi*Re))
        berg.X += dt*berg.Vx/(np.cos((((y_ + berg.Y)/2)*np.pi)/180))*(180/(np.pi*Re))
        
        # update date iceberg history if still in bounds
        if berg.in_bounds(x_bounds, y_bounds):
            berg.update_history()
        
        # otherwise, terminate
        else:
            berg.out_of_bounds = True
            break
      
    return berg

def ab2(berg, ocean_data, atm_data, drift, t_all):
    
    Re = 6378*1e3  # radius of Earth
    
    # determine bounding box for drift simulation
    x_bounds, y_bounds = get_bounding_box(ocean_data, atm_data)
    
    
    # get timestep
    tdelta = t_all[1] - t_all[0]
    dt = tdelta.seconds
    
    # drift iceberg for all t in t_all
    for i, t in enumerate(t_all):
        
        berg.T = t
    
        Vcx, Vcy, Vax, Vay = get_interpolated_values(berg, ocean_data, atm_data)
            
        # integrate
        if i < 1:
            # explicit Euler forward scheme
            berg.Vx += dt*berg.Ax
            berg.Vy += dt*berg.Ay   
            berg.Ax, berg.Ay = drift(berg, Vax, Vay, Vcx, Vcy)
            
        else:
            # second order Adams Bashforth
            berg.Vx += dt*(1.5*berg.Ax - 0.5*berg.history['Ax'][-1])
            berg.Vy += dt*(1.5*berg.Ay - 0.5*berg.history['Ay'][-1])
            berg.Ax, berg.Ay = drift(berg, Vax, Vay, Vcx, Vcy)
         
        # calculate position in degrees lat/lon
        y_ = berg.Y
        berg.Y += dt*berg.Vy*(180/(np.pi*Re))
        berg.X += dt*berg.Vx/(np.cos((((y_ + berg.Y)/2)*np.pi)/180))*(180/(np.pi*Re))
        
        # update date iceberg history if still in bounds
        if berg.in_bounds(x_bounds, y_bounds):
            berg.update_history()
        
        # otherwise, terminate
        else:
            berg.out_of_bounds = True
            break
      
    
    return berg

def ab3(berg, ocean_data, atm_data, drift, t_all):
    
    Re = 6378*1e3  # radius of Earth
    
    # determine bounding box for drift simulation
    x_bounds, y_bounds = get_bounding_box(ocean_data, atm_data)
    
    # get timestep
    tdelta = t_all[1] - t_all[0]
    dt = tdelta.seconds
    
    # drift iceberg for all t in t_all
    for i, t in enumerate(t_all):
        
        berg.T = t
    
        Vcx, Vcy, Vax, Vay = get_interpolated_values(berg, ocean_data, atm_data)
            
        # integrate
        if i < 1:
            # explicit Euler forward scheme
            berg.Vx += dt*berg.Ax
            berg.Vy += dt*berg.Ay   
            berg.Ax, berg.Ay = drift(berg, Vax, Vay, Vcx, Vcy)
            
        elif i < 3:
            # second order Adams Bashforth
            berg.Vx += dt*(1.5*berg.Ax - 0.5*berg.history['Ax'][-1])
            berg.Vy += dt*(1.5*berg.Ay - 0.5*berg.history['Ay'][-1])
            berg.Ax, berg.Ay = drift(berg, Vax, Vay, Vcx, Vcy)
         
        else:
            # third order Adams Bashforth
            berg.Vx += dt/12*(23*berg.Ax-16*berg.history['Ax'][-1]+5*berg.history['Ax'][-2])
            berg.Vy += dt/12*(23*berg.Ay-16*berg.history['Ay'][-1]+5*berg.history['Ay'][-2])
            berg.Ax, berg.Ay = drift(berg, Vax, Vay, Vcx, Vcy)
            
            
        # calculate position in degrees lat/lon
        y_ = berg.Y
        berg.Y += dt*berg.Vy*(180/(np.pi*Re))
        berg.X += dt*berg.Vx/(np.cos((((y_ + berg.Y)/2)*np.pi)/180))*(180/(np.pi*Re))
        
        # update date iceberg history if still in bounds
        if berg.in_bounds(x_bounds, y_bounds):
            berg.update_history()
        
        # otherwise, terminate
        else:
            berg.out_of_bounds = True
            break
      
    
    return berg

def ab4(berg, ocean_data, atm_data, drift, t_all):
    
    Re = 6378*1e3  # radius of Earth
    
    # determine bounding box for drift simulation
    x_bounds, y_bounds = get_bounding_box(ocean_data, atm_data)
    
    
    # get timestep
    tdelta = t_all[1] - t_all[0]
    dt = tdelta.seconds
    
    # drift iceberg for all t in t_all
    for i, t in enumerate(t_all):
        
        berg.T = t
    
        Vcx, Vcy, Vax, Vay = get_interpolated_values(berg, ocean_data, atm_data)
            
        # integrate
        if i < 1:
            # explicit Euler forward scheme
            berg.Vx += dt*berg.Ax
            berg.Vy += dt*berg.Ay   
            berg.Ax, berg.Ay = drift(berg, Vax, Vay, Vcx, Vcy)
            
        elif i < 3:
            # second order Adams-Bashforth
            berg.Vx += dt*(1.5*berg.Ax - 0.5*berg.history['Ax'][-1])
            berg.Vy += dt*(1.5*berg.Ay - 0.5*berg.history['Ay'][-1])
            berg.Ax, berg.Ay = drift(berg, Vax, Vay, Vcx, Vcy)
         
        else:
            # fourth order Adams-Bashforth predictor-corrector
            berg.Vx += dt/24*(55*berg.Ax-59*berg.history['Ax'][-1]+37*berg.history['Ax'][-2]-9*berg.history['Ax'][-3])
            berg.Vy += dt/24*(55*berg.Ay-59*berg.history['Ay'][-1]+37*berg.history['Ay'][-2]-9*berg.history['Ay'][-3])
            
            Ax1, Ay1 = drift(berg, Vax, Vay, Vcx, Vcy)
            
            berg.Vx += dt/24*(Ax1+19*berg.Ax-5*berg.history['Ax'][-1]+berg.history['Ax'][-2])
            berg.Vy += dt/24*(Ay1+19*berg.Ay-5*berg.history['Ay'][-1]+berg.history['Ay'][-2])
            
            berg.Ax, berg.Ay = drift(berg, Vax, Vay, Vcx, Vcy)

            
        # calculate position in degrees lat/lon
        y_ = berg.Y
        berg.Y += dt*berg.Vy*(180/(np.pi*Re))
        berg.X += dt*berg.Vx/(np.cos((((y_ + berg.Y)/2)*np.pi)/180))*(180/(np.pi*Re))
        
        # update date iceberg history if still in bounds
        if berg.in_bounds(x_bounds, y_bounds):
            berg.update_history()
        
        # otherwise, terminate
        else:
            berg.out_of_bounds = True
            break
      
    
    return berg