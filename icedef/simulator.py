"""This module does stuff.

"""

#import numpy as np
#import xarray as xr
#from scipy.optimize import minimize
#from icedef import iceberg, metocean, drift, tools


def run_optimization(reference_vectors, start_location, time_frame):

    optimization_result = minimize(optimization_wrapper, x0=(1, 1), bounds=((0, 15), (0, 15)),
                                   args=(reference_vectors, start_location, time_frame))

    return optimization_result


def compute_mse(simulation_vectors, reference_vectors):

    sim_lats, sim_lons = simulation_vectors
    ref_lats, ref_lons = reference_vectors

    mean_square_error_list = []

    for i in range(len(ref_lats)):

        time = ref_lats['time'][i]
        sim_lat = sim_lats.interp(time=time)
        sim_lon = sim_lons.interp(time=time)
        mean_square_error_list.append(np.sqrt((sim_lat - ref_lats[i])**2 + (sim_lon - ref_lons[i])**2))

    return np.mean(mean_square_error_list)


def optimization_wrapper(drag_coeffs, reference_vectors, start_location, time_frame):

    Ca, Cw = drag_coeffs
    xds = run_simulation(start_location, time_frame, Ca=Ca, Cw=Cw)
    simulation_vectors = (xds['latitude'], xds['longitude'])
    mse = compute_mse(simulation_vectors, reference_vectors)

    return mse


def run_simulation(start_location, time_frame, **kwargs):

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

    times = np.zeros(nt, dtype='datetime64[ns]')

    results = {'latitude': np.zeros(nt),
               'longitude': np.zeros(nt),
               'iceberg_eastward_velocity': np.zeros(nt),
               'iceberg_northward_velocity': np.zeros(nt),
               'current_eastward_velocity': np.zeros(nt),
               'current_northward_velocity': np.zeros(nt),
               'wind_eastward_velocity': np.zeros(nt),
               'wind_northward_velocity': np.zeros(nt),
               'current_eastward_force': np.zeros(nt),
               'current_northward_force': np.zeros(nt),
               'wind_eastward_force': np.zeros(nt),
               'wind_northward_force': np.zeros(nt),
               'coriolis_eastward_force': np.zeros(nt),
               'coriolis_northward_force': np.zeros(nt),
               'pressure_eastward_force': np.zeros(nt),
               'pressure_northward_force': np.zeros(nt)}

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
            'latitude': iceberg_.latitude,
            'ekman': kwargs.pop('ekman', False),
            'depth_vec': kwargs.pop('depth_vec', np.arange(0, -110, -10))
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
                                                       **iceberg_constants)

                # Store results from this timestep

                times[i] = iceberg_.time
                results['latitude'][i] = iceberg_.latitude
                results['longitude'][i] = iceberg_.longitude
                results['iceberg_eastward_velocity'][i] = iceberg_.eastward_velocity
                results['iceberg_northward_velocity'][i] = iceberg_.northward_velocity
                results['current_eastward_velocity'][i] = current_velocity[0]
                results['current_northward_velocity'][i] = current_velocity[1]
                results['wind_eastward_velocity'][i] = wind_velocity[0]
                results['wind_northward_velocity'][i] = wind_velocity[1]
                results['wind_eastward_force'][i] = forces[0]
                results['wind_northward_force'][i] = forces[1]
                results['current_eastward_force'][i] = forces[2]
                results['current_northward_force'][i] = forces[3]
                results['coriolis_eastward_force'][i] = forces[4]
                results['coriolis_northward_force'][i] = forces[5]
                results['pressure_eastward_force'][i] = forces[6]
                results['pressure_northward_force'][i] = forces[7]

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

    xds = xr.Dataset()

    for key, value in results.items():

        xarr = xr.DataArray(data=value, coords=[times], dims=['time'])
        xds[key] = xarr

    return xds
