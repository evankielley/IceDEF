"""This module sets up and runs iceberg drift simulations and optimizations.
"""

import numpy as np
import xarray as xr
from copy import deepcopy
from scipy.optimize import minimize
from icedef import iceberg, metocean, drift, tools, timesteppers


class Simulator:

    def __init__(self, **kwargs):
        """This class sets up and runs the components necessary to run an iceberg drift simulation.

        Kwargs:
            time_step (numpy.timedelta64): Time step in seconds.
            drift_model (function): The drift model function.
            time_stepper (function): The numerical integrator function.
            ocean_model (str): Name of ocean model. Can be ECMWF or HYCOM.
            atmosphere_model (str): Name of the atmosphere model. Can be ECMWF or NARR.
        """

        self.time_step = kwargs.pop('time_step', np.timedelta64(300, 's'))
        self.drift_model = kwargs.pop('drift_model', drift.newtonian_drift_wrapper)
        self.time_stepper = kwargs.pop('time_stepper', timesteppers.euler)
        self.ocean_model = kwargs.pop('ocean_model', 'ECMWF')
        self.atmosphere_model = kwargs.pop('atmosphere_model', 'NARR')
        self.start_latitude = kwargs.pop('start_latitude', None)
        self.start_longitude = kwargs.pop('start_longitude', None)
        self._start_location = None
        self.start_time = kwargs.pop('start_time', None)
        self.end_time = kwargs.pop('end_time', None)
        self._time_frame = None
        self.results = {}

    @property
    def start_location(self):
        return self.start_latitude, self.start_longitude

    @property
    def time_frame(self):
        return self.start_time, self.end_time

    def run_simulation(self, store_results_as=None, **kwargs):
        """This method simulates iceberg drift.

        Args:
            store_results_as (str): Key by which the results of the simulation will be saved in results attribute.
        """

        kwargs['time_step'] = self.time_step
        kwargs['time_stepper'] = self.time_stepper
        kwargs['drift_model'] = self.drift_model
        kwargs['ocean_model'] = self.ocean_model
        kwargs['atmosphere_model'] = self.atmosphere_model

        results = run_simulation(self.start_location, self.time_frame, **kwargs)

        if store_results_as is not None:

            self.results[store_results_as] = results

        return results

    def run_optimization(self, keys, x0, bounds, reference_vectors):
        """This function optimizes user specified drift simulation parameters using the Scipy minimize function.

        Args:
            keys (list of str): The names of the drift simulation kwargs to be optimized.
            x0 (numpy.ndarray): The initial guesses.
            bounds (list of list of float): The upper and lower bounds for the parameters being optimized.
            reference_vectors (tuple of xarray.core.dataarray.DataArray): The latitude, longitude vectors to compare to.

        Returns:
            optimization_result (scipy.optimize.optimize.OptimizeResult): Results from minimization.

        """

        optimization_result = minimize(self.optimization_wrapper, x0=x0, bounds=bounds, args=(keys, reference_vectors))

        return optimization_result

    def optimization_wrapper(self, values, keys, reference_vectors):

        kwargs = dict(zip(keys, values))
        xds = self.run_simulation(**kwargs)
        simulation_vectors = (xds['latitude'], xds['longitude'])
        mean_square_errors = compute_mse(simulation_vectors, reference_vectors)

        return np.mean(mean_square_errors)


def compute_mse(simulation_vectors, reference_vectors):

    sim_lats, sim_lons = simulation_vectors
    ref_lats, ref_lons = reference_vectors

    mean_square_errors = []

    stop_index = np.where(ref_lats['time'].values <= sim_lats['time'].values[-1])[0][-1]

    for i in range(stop_index + 1):

        time = ref_lats['time'][i]
        sim_lat = sim_lats.interp(time=time, assume_sorted=True)
        sim_lon = sim_lons.interp(time=time, assume_sorted=True)
        mean_square_error = np.sqrt((sim_lat - ref_lats[i])**2 + (sim_lon - ref_lons[i])**2)
        mean_square_errors.append(mean_square_error)

    return mean_square_errors


def run_simulation(start_location, time_frame, **kwargs):

    time_step = kwargs.pop('time_step', np.timedelta64(300, 's'))
    time_stepper = kwargs.pop('time_stepper', timesteppers.euler)
    drift_model = kwargs.pop('drift_model', drift.newtonian_drift_wrapper)
    ocean_model = kwargs.pop('ocean_model', 'ECMWF')
    atmosphere_model = kwargs.pop('atmosphere_model', 'NARR')

    start_time, end_time = time_frame
    dt = time_step.item().total_seconds()
    nt = int((end_time - start_time).item().total_seconds() / dt)

    waterline_length = kwargs.pop('waterline_length', 160)
    sail_height = kwargs.pop('sail_height', 60)
    size = waterline_length, sail_height
    start_velocity = kwargs.pop('start_velocity', (0, 0))

    iceberg_ = iceberg.quickstart(start_time, start_location, velocity=start_velocity, size=size, **kwargs)

    ocean = metocean.Ocean(time_frame, model=ocean_model)
    atmosphere = metocean.Atmosphere(time_frame, model=atmosphere_model)

    # Initialize arrays
    times = np.zeros(nt, dtype='datetime64[ns]')
    results = {'latitude': np.zeros(nt),
               'longitude': np.zeros(nt),
               'easting': np.zeros(nt),
               'northing': np.zeros(nt),
               'iceberg_eastward_velocity': np.zeros(nt),
               'iceberg_northward_velocity': np.zeros(nt)}
    kwargs = {
        'form_drag_coefficient_in_air': kwargs.pop('Ca', iceberg_.FORM_DRAG_COEFFICIENT_IN_AIR),
        'form_drag_coefficient_in_water': kwargs.pop('Cw', iceberg_.FORM_DRAG_COEFFICIENT_IN_WATER),
        'skin_drag_coefficient_in_air': iceberg_.SKIN_DRAG_COEFFICIENT_IN_AIR,
        'skin_drag_coefficient_in_water': iceberg_.SKIN_DRAG_COEFFICIENT_IN_WATER,
        'sail_area': iceberg_.geometry.sail_area,
        'keel_area': iceberg_.geometry.keel_area,
        'top_area': iceberg_.geometry.waterline_length ** 2,
        'bottom_area': iceberg_.geometry.bottom_area,
        'mass': kwargs.pop('mass', iceberg_.geometry.mass),
        'latitude': iceberg_.latitude,
        'ekman': kwargs.pop('ekman', False),
        'depth_vec': kwargs.pop('depth_vec', np.arange(0, -110, -10)),
        'current_acceleration': (0, 0),
        'eastward_current': ocean.current.eastward_velocities,
        'northward_current': ocean.current.northward_velocities,
        'eastward_wind': atmosphere.wind.eastward_velocities,
        'northward_wind': atmosphere.wind.northward_velocities,
    }

    for i in range(nt):

        times[i] = iceberg_.time
        results['latitude'][i] = iceberg_.latitude
        results['longitude'][i] = iceberg_.longitude
        results['easting'][i] = iceberg_.easting
        results['northing'][i] = iceberg_.northing
        results['iceberg_eastward_velocity'][i] = iceberg_.eastward_velocity
        results['iceberg_northward_velocity'][i] = iceberg_.northward_velocity

        if time_stepper in (timesteppers.ab2, timesteppers.ab3):

            dx, dy, dvx, dvy = time_stepper(drift_model, dt,
                                            times[:i+1],
                                            results['longitude'][:i+1],
                                            results['latitude'][:i+1],
                                            results['iceberg_eastward_velocity'][:i+1],
                                            results['iceberg_northward_velocity'][:i+1],
                                            **kwargs)

        else:

            dx, dy, dvx, dvy = time_stepper(drift_model, dt,
                                            iceberg_.time, iceberg_.longitude, iceberg_.latitude,
                                            iceberg_.eastward_velocity, iceberg_.northward_velocity,
                                            **kwargs)

        iceberg_.eastward_velocity += dvx
        iceberg_.northward_velocity += dvy
        iceberg_.easting += dx
        iceberg_.northing += dy
        iceberg_.time += time_step
        iceberg_.latitude += tools.dy_to_dlat(dy)
        iceberg_.longitude += tools.dx_to_dlon(dx, iceberg_.latitude)

    xds = xr.Dataset()

    for key, value in results.items():
        xarr = xr.DataArray(data=value, coords=[times], dims=['time'])
        xds[key] = xarr

    return xds


def run_test_simulation(start_location, time_frame, **kwargs):
    # Kwargs for testing
    constant_current_velocity = kwargs.pop('constant_current_velocity', None)
    constant_eastward_current_velocity = kwargs.pop('constant_eastward_current_velocity', None)
    constant_northward_current_velocity = kwargs.pop('constant_northward_current_velocity', None)
    constant_wind_velocity = kwargs.pop('constant_wind_velocity', None)
    constant_eastward_wind_velocity = kwargs.pop('constant_eastward_wind_velocity', None)
    constant_northward_wind_velocity = kwargs.pop('constant_northward_wind_velocity', None)
    constant_current_acceleration = kwargs.pop('constant_current_acceleration', False)

    # Args
    start_latitude, start_longitude = start_location
    start_time, end_time = time_frame

    # Kwargs
    start_velocity = kwargs.pop('start_velocity', (0, 0))
    time_step = kwargs.pop('time_step', np.timedelta64(300, 's'))
    drift_model = kwargs.pop('drift_model', 'Newtonian')
    numerical_method = kwargs.pop('numerical_method', 'Euler')

    # Object creation
    iceberg_ = iceberg.quickstart(start_time, start_location, velocity=start_velocity, **kwargs)
    metocean_ = metocean.Metocean(time_frame)
    ocean = metocean_.ocean
    atmosphere = metocean_.atmosphere

    assume_constant_current_velocity = False

    if constant_eastward_current_velocity is not None and constant_northward_current_velocity is not None:
        assume_constant_current_velocity = True
        constant_current_velocity = (constant_eastward_current_velocity, constant_northward_current_velocity)
        current_velocity = constant_current_velocity

    if constant_current_velocity is not None:
        assume_constant_current_velocity = True
        current_velocity = constant_current_velocity

    elif constant_eastward_current_velocity is not None:
        ocean.eastward_velocities.values = np.full(ocean.eastward_velocities.shape,
                                                           constant_eastward_current_velocity)

    elif constant_northward_current_velocity is not None:
        ocean.northward_current_velocities.values = np.full(ocean.northward_current_velocities.shape,
                                                            constant_northward_current_velocity)

    if not assume_constant_current_velocity:
        current_velocity_interpolator = metocean.Interpolate((ocean.data.time.values,
                                                              ocean.data.latitude.values,
                                                              ocean.data.longitude.values),
                                                             ocean.eastward_velocities.values,
                                                             ocean.northward_current_velocities.values)

    assume_constant_wind_velocity = False

    if constant_eastward_wind_velocity is not None and constant_northward_wind_velocity is not None:
        assume_constant_wind_velocity = True
        constant_wind_velocity = (constant_eastward_wind_velocity, constant_northward_wind_velocity)
        wind_velocity = constant_wind_velocity

    if constant_wind_velocity is not None:
        assume_constant_wind_velocity = True
        wind_velocity = constant_wind_velocity

    elif constant_eastward_wind_velocity is not None:
        atmosphere.eastward_wind_velocities.values = np.full(atmosphere.eastward_wind_velocities.shape,
                                                             constant_eastward_wind_velocity)

    elif constant_northward_wind_velocity is not None:
        atmosphere.northward_wind_velocities.values = np.full(atmosphere.northward_wind_velocities.shape,
                                                              constant_northward_wind_velocity)

    if not assume_constant_wind_velocity:
        wind_velocity_interpolator = metocean.Interpolate((atmosphere.data.time.values,
                                                           atmosphere.data.latitude.values,
                                                           atmosphere.data.longitude.values),
                                                          atmosphere.eastward_wind_velocities.values,
                                                          atmosphere.northward_wind_velocities.values)

    dt = time_step.item().total_seconds()
    nt = int((end_time - start_time).item().total_seconds() / dt)

    times = np.zeros(nt, dtype='datetime64[ns]')

    results = {'latitude': np.zeros(nt),
               'longitude': np.zeros(nt),
               'easting': np.zeros(nt),
               'northing': np.zeros(nt),
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

        constants = {
            'form_drag_coefficient_in_air': kwargs.pop('Ca', iceberg_.FORM_DRAG_COEFFICIENT_IN_AIR),
            'form_drag_coefficient_in_water': kwargs.pop('Cw', iceberg_.FORM_DRAG_COEFFICIENT_IN_WATER),
            'skin_drag_coefficient_in_air': iceberg_.SKIN_DRAG_COEFFICIENT_IN_AIR,
            'skin_drag_coefficient_in_water': iceberg_.SKIN_DRAG_COEFFICIENT_IN_WATER,
            'sail_area': iceberg_.geometry.sail_area,
            'keel_area': iceberg_.geometry.keel_area,
            'top_area': iceberg_.geometry.waterline_length ** 2,
            'bottom_area': 0,
            'mass': iceberg_.geometry.mass,
            'latitude': iceberg_.latitude,
            'ekman': kwargs.pop('ekman', False),
            'depth_vec': kwargs.pop('depth_vec', np.arange(0, -110, -10)),
            'current_acceleration': (0, 0),
            'zero_current_force': kwargs.pop('zero_current_force', False),
            'zero_wind_force': kwargs.pop('zero_wind_force', False),
            'zero_coriolis_force': kwargs.pop('zero_coriolis_force', False),
            'zero_pressure_force': kwargs.pop('zero_pressure_force', False),
        }

        point = (iceberg_.time, iceberg_.latitude, iceberg_.longitude)

        if not assume_constant_current_velocity:
            current_velocity = current_velocity_interpolator.interpolate(point)

        if not assume_constant_wind_velocity:
            wind_velocity = wind_velocity_interpolator.interpolate(point)

        if numerical_method == 'Euler':

            for i in range(nt):

                # Compute instantaneous acceleration from drift model
                ax, ay, forces = drift.newtonian_drift((iceberg_.eastward_velocity, iceberg_.northward_velocity),
                                                       current_velocity, wind_velocity, **constants)

                # Store results from this time step
                times[i] = iceberg_.time
                results['latitude'][i] = iceberg_.latitude
                results['longitude'][i] = iceberg_.longitude
                results['easting'][i] = iceberg_.easting
                results['northing'][i] = iceberg_.northing
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

                # Update parameters for next time step
                iceberg_.time += time_step
                iceberg_.eastward_velocity += ax * dt
                iceberg_.northward_velocity += ay * dt
                dy = iceberg_.northward_velocity * dt
                iceberg_.northing += dy
                iceberg_.latitude += tools.dy_to_dlat(dy)
                constants['latitude'] = iceberg_.latitude
                dx = iceberg_.eastward_velocity * dt
                iceberg_.easting += dx
                iceberg_.longitude += tools.dx_to_dlon(dx, iceberg_.latitude)

                point = (iceberg_.time, iceberg_.latitude, iceberg_.longitude)

                old_current_velocity = deepcopy(current_velocity)

                if not assume_constant_current_velocity:
                    current_velocity = current_velocity_interpolator.interpolate(point)

                if not assume_constant_wind_velocity:
                    wind_velocity = wind_velocity_interpolator.interpolate(point)

                if constant_current_acceleration is not False:
                    constants['current_acceleration'] = constant_current_acceleration

                else:
                    current_acceleration = ((current_velocity[0] - old_current_velocity[0]) / dt,
                                            (current_velocity[1] - old_current_velocity[1]) / dt)
                    constants['current_acceleration'] = current_acceleration

    xds = xr.Dataset()

    for key, value in results.items():
        xarr = xr.DataArray(data=value, coords=[times], dims=['time'])
        xds[key] = xarr

    return xds
