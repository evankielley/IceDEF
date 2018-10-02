import numpy as np
from icedef import iceberg, metocean, drift, tools

class Simulator:

    def __init__(self, **kwargs):
        
        self.iceberg = kwargs.pop('iceberg', None)
        self.metocean = kwargs.pop('metocean', None)
        self.results = None
    
    
    def run_simulation(self, start_location, time_frame, **kwargs):
        
        # Args
        start_latitude, start_longitude = start_location
        start_time, end_time = time_frame
        
        # Kwargs
        start_velocity = kwargs.pop('start_velocity', (0, 0))
        time_step = kwargs.pop('time_step', np.timedelta64(300, 's'))
        drift_model = kwargs.pop('drift_model', 'Newtonian')
        numerical_method = kwargs.pop('numerical_method', 'Euler')
        
        # Object creation
        iceberg_ = iceberg.quickstart(start_time, start_location, velocity=start_velocity)
        metocean_ = metocean.Metocean(time_frame)
        
        
        dt = time_step.item().total_seconds()
        nt = int((end_time - start_time).item().total_seconds() / dt)
        
        results = {'time': [None] * nt,
                   'iceberg_position': [None] * nt, 
                   'iceberg_velocity': [None] * nt,
                   'current_velocity': [None] * nt,
                   'wind_velocity': [None] * nt,
                   'current_force': [None] * nt,
                   'wind_force': [None] * nt,
                   'coriolis_force': [None] * nt,
                   'pressure_gradient_force': [None] * nt}
        
        if drift_model == 'Newtonian':
            
            iceberg_constants = {
                'form_drag_coefficient_in_air': iceberg_.FORM_DRAG_COEFFICIENT_IN_AIR,
                'form_drag_coefficient_in_water': iceberg_.FORM_DRAG_COEFFICIENT_IN_WATER,
                'skin_drag_coefficient_in_air': iceberg_.SKIN_DRAG_COEFFICIENT_IN_AIR,
                'skin_drag_coefficient_in_water': iceberg_.SKIN_DRAG_COEFFICIENT_IN_WATER,
                'sail_area': iceberg_.geometry.sail_area,
                'keel_area': iceberg_.geometry.keel_area,
                'top_area': iceberg_.geometry.waterline_length**2,
                'bottom_area': 0,
                'mass': iceberg_.geometry.mass,
                'latitude': iceberg_.latitude
                }
            
            point = (iceberg_.time, iceberg_.latitude, iceberg_.longitude)
    
            current_velocity = (metocean_.interpolate(point, metocean_.ocean.eastward_current_velocities),
                                metocean_.interpolate(point, metocean_.ocean.northward_current_velocities))

            wind_velocity = (metocean_.interpolate(point, metocean_.atmosphere.eastward_wind_velocities),
                             metocean_.interpolate(point, metocean_.atmosphere.northward_wind_velocities))
            
            
            if numerical_method == 'Euler':
                                
                for i in range(nt):
                    
                    # Compute instantaneous acceleration from drift model    
                    ax, ay, forces = drift.newtonian_drift((iceberg_.eastward_velocity, iceberg_.northward_velocity), 
                                                           current_velocity, wind_velocity, 
                                                           iceberg_constants)

                    # Store results from this timestep
                    
                    results['time'][i] = iceberg_.time
                    results['iceberg_position'][i] = iceberg_.latitude, iceberg_.longitude
                    results['iceberg_velocity'][i] = iceberg_.eastward_velocity, iceberg_.northward_velocity
                    results['current_velocity'][i] = current_velocity
                    results['wind_velocity'][i] = wind_velocity
                    results['current_force'][i] = forces[2], forces[3]
                    results['wind_force'][i] = forces[0], forces[1]
                    results['coriolis_force'][i] = forces[4], forces[5]
                    results['pressure_gradient_force'][i] = forces[6], forces[7]
                    
                    # Update parameters for next timestep
                    
                    iceberg_.time += time_step
                    iceberg_.eastward_velocity += ax * dt
                    iceberg_.northward_velocity += ay * dt
                    iceberg_.latitude += tools.dy_to_dlat(iceberg_.northward_velocity * dt)
                    iceberg_constants['latitude'] = iceberg_.latitude
                    iceberg_.longitude += tools.dx_to_dlon(iceberg_.eastward_velocity * dt, iceberg_.latitude)
                    point = (iceberg_.time, iceberg_.latitude, iceberg_.longitude)
                    current_velocity = (metocean_.interpolate(point, metocean_.ocean.eastward_current_velocities),
                                        metocean_.interpolate(point, metocean_.ocean.northward_current_velocities))
                    wind_velocity = (metocean_.interpolate(point, metocean_.atmosphere.eastward_wind_velocities),
                                     metocean_.interpolate(point, metocean_.atmosphere.northward_wind_velocities))

                                                                    
        self.results = results
                      
        
        