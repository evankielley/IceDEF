import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from icedef import iceberg, metocean, drift, tools

class Simulator:

    def __init__(self, **kwargs):
        
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
        
        results = {'time': np.zeros(nt, dtype='datetime64[ns]'),
                   'iceberg_position': np.zeros((nt, 2)), 
                   'iceberg_velocity': np.zeros((nt, 2)),
                   'current_velocity': np.zeros((nt, 2)),
                   'wind_velocity': np.zeros((nt, 2)),
                   'current_force': np.zeros((nt, 2)),
                   'wind_force': np.zeros((nt, 2)),
                   'coriolis_force': np.zeros((nt, 2)),
                   'pressure_gradient_force': np.zeros((nt, 2))}
        
        if drift_model == 'Newtonian':
            
            iceberg_constants = {
                'form_drag_coefficient_in_air': kwargs.pop('Ca', iceberg_.FORM_DRAG_COEFFICIENT_IN_AIR),
                'form_drag_coefficient_in_water': kwargs.pop('Cw', iceberg_.FORM_DRAG_COEFFICIENT_IN_WATER),
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

                                                                    
        return results
        
        
    def run_optimization(self, reference_points, start_location, time_frame):
        
        optimization_result = minimize(self.optimization_wrapper, x0=(1, 1), bounds=((0, 15), (0, 15)),
                                       args=(reference_points, start_location, time_frame))
        
        return optimization_result
            
        
    def optimization_wrapper(self, drag_coeffs, reference_points, start_location, time_frame):
        
        Ca, Cw = drag_coeffs
        results = self.run_simulation(start_location, time_frame, Ca=Ca, Cw=Cw)
        mse = self.compute_mse((results['iceberg_position'][:, 1], results['iceberg_position'][:, 0], results['time']),
                               reference_points)
        return mse
            
    def compute_mse(self, simulation_vectors, reference_points):
        
        t0 = np.datetime64('1900-01-01T00:00:00')
        
        sim_x_vec, sim_y_vec, sim_t_vec = simulation_vectors
        #print(type((np.datetime64(sim_t_vec[0]) - t0)))
        #print(type((np.datetime64(sim_t_vec[0]) - t0).item()))
        #print(type((np.datetime64(sim_t_vec[0]) - t0).total_seconds()))
        sim_t_vec = [(np.datetime64(sim_t) - t0).item() for sim_t in sim_t_vec]

        ref_x_vec, ref_y_vec, ref_t_vec = reference_points
        
        mse_list = []
        
        for i in range(len(ref_t_vec)):
            
            ref_t = ref_t_vec[i]
            ref_x = ref_x_vec[i]
            ref_y = ref_y_vec[i]
            
            ref_t = (ref_t - t0).item()

            sim_x_interpolant = interp1d(sim_t_vec, sim_x_vec)
            sim_y_interpolant = interp1d(sim_t_vec, sim_y_vec)

            sim_x = sim_x_interpolant(ref_t)
            sim_y = sim_y_interpolant(ref_t)
            
            mse = np.sqrt((ref_x - sim_x)**2 + (ref_y - sim_y)**2)
            mse_list.append(mse)
            
        mean_mse = np.mean(np.array(mse_list))
        
        return mean_mse
        
        
                      
    
        