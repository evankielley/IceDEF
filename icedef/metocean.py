import os
import numpy as np
import xarray as xr
from urllib.request import urlretrieve
from datetime import date, timedelta
from pandas import Timestamp


class Metocean:

    def __init__(self, date_bounds, **kwargs):

        self.date_bounds = date_bounds
        self.ocean_model = kwargs.pop('ocean_model', 'ECMWF')
        self.atmosphere_model = kwargs.pop('atmosphere_model', 'NARR')
        self.ocean = Ocean(date_bounds, model=self.ocean_model)
        self.atmosphere = Atmosphere(date_bounds, model=self.atmosphere_model)

    def swap_ocean(self, ocean_model):

        self.ocean_model = ocean_model
        self.ocean = self.Ocean(self.date_bounds, model=ocean_model)

    def swap_atmosphere(self, atmosphere_model):

        self.atmosphere_model = atmosphere_model
        self.atmosphere = self.Atmosphere(self.date_bounds, model=atmosphere_model)


class Ocean:

    ID = None
    PATH = None

    def __init__(self, date_bounds, model='ECMWF'):

        if model == 'ECMWF':

            self.ID = "GLOBAL_ANALYSIS_FORECAST_PHY_001_024"
            self.PATH = 'http://icedef.munroelab.ca/data/ECMWF/ocean/daily/'
            self.data = xr.open_mfdataset(get_files(self.ID, self.PATH, date_bounds)).squeeze('depth').rename(
                                                    {'uo': 'eastward_velocity', 'vo': 'northward_velocity'})

        elif model == 'HYCOM':

            self.ID = "HYCOM_Region_1_3D"
            self.PATH = 'http://icedef.munroelab.ca/data/HYCOM/ocean/daily/'
            self.data = xr.open_mfdataset(get_files(self.ID, self.PATH, date_bounds)).rename(
                {'eastward_current_velocity': 'eastward_velocity', 'northward_current_velocity': 'northward_velocity'})

        else:

            print('Invalid model.')

        self.current = self.Current(self.data)

    class Current:

        def __init__(self, data):

            self.eastward_velocities = xr.DataArray(data=data.eastward_velocity.values,
                                                    coords=[('time', data.time.values),
                                                            ('latitude', data.latitude.values),
                                                            ('longitude', data.longitude.values)],
                                                    attrs=data.eastward_velocity.attrs)

            self.northward_velocities = xr.DataArray(data=data.northward_velocity.values,
                                                     coords=[('time', data.time.values),
                                                             ('latitude', data.latitude.values),
                                                             ('longitude', data.longitude.values)],
                                                     attrs=data.eastward_velocity.attrs)

            self.speeds = xr.DataArray(data=compute_magnitude(self.eastward_velocities,
                                                              self.northward_velocities),
                                       coords=[('time', data.time.values),
                                               ('latitude', data.latitude.values),
                                               ('longitude', data.longitude.values)])

            self.directions = xr.DataArray(data=compute_direction(self.eastward_velocities,
                                                                  self.northward_velocities),
                                           coords=[('time', data.time.values),
                                                   ('latitude', data.latitude.values),
                                                   ('longitude', data.longitude.values)])


class Atmosphere:

    ID = None
    PATH = None

    def __init__(self, date_bounds, model='NARR'):

        if model == 'ECMWF':

            self.ID = "WIND_GLO_WIND_L4_NRT_OBSERVATIONS_012_004"
            self.PATH = 'http://icedef.munroelab.ca/data/ECMWF/atm/daily/'
            self.data = xr.open_mfdataset(get_files(self.ID, self.PATH, date_bounds)).squeeze('depth').rename(
                {'eastward_wind': 'eastward_velocity',
                 'northward_wind': 'northward_velocity'})

        elif model == 'NARR':

            self.ID = "NCEP_North_American_Regional_Reanalysis_NARR"
            self.PATH = 'http://icedef.munroelab.ca/data/NARR/atm/daily/'
            self.data = xr.open_mfdataset(get_files(self.ID, self.PATH, date_bounds)).rename(
                {'uwnd': 'eastward_velocity',
                 'vwnd': 'northward_velocity',
                 'lat': 'latitude',
                 'lon': 'longitude'})
            self.data['longitude'] = np.mod(self.data.longitude - 180, 360) - 180

        else:

            print('Invalid model.')

        self.wind = self.Wind(self.data)

    class Wind:

        def __init__(self, data):

            self.eastward_velocities = xr.DataArray(data=data.eastward_velocity.values,
                                                    coords=[('time', data.time.values),
                                                            ('latitude', data.latitude.values),
                                                            ('longitude', data.longitude.values)],
                                                    attrs=data.eastward_velocity.attrs)

            self.northward_velocities = xr.DataArray(data=data.northward_velocity.values,
                                                     coords=[('time', data.time.values),
                                                             ('latitude', data.latitude.values),
                                                             ('longitude', data.longitude.values)],
                                                     attrs=data.eastward_velocity.attrs)

            self.speeds = xr.DataArray(data=compute_magnitude(self.eastward_velocities,
                                                              self.northward_velocities),
                                       coords=[('time', data.time.values),
                                               ('latitude', data.latitude.values),
                                               ('longitude', data.longitude.values)])

            self.directions = xr.DataArray(data=compute_direction(self.eastward_velocities,
                                                                  self.northward_velocities),
                                           coords=[('time', data.time.values),
                                                   ('latitude', data.latitude.values),
                                                   ('longitude', data.longitude.values)])


def get_files(id_, path, date_bounds, cache=True):

    start_date, end_date = date_bounds

    if isinstance(start_date, np.datetime64):
        start_date = Timestamp(start_date)

    if isinstance(end_date, np.datetime64):
        end_date = Timestamp(end_date)

    start_date = date(start_date.year, start_date.month, start_date.day)
    end_date = date(end_date.year, end_date.month, end_date.day)
    time_delta = end_date - start_date

    if cache:
        cache_path = 'cache/' + id_ + '/'
        if not os.path.exists(cache_path):
            try:
                os.makedirs(cache_path)
            except:
                print("Couldn't make directory for cache")
    else:
        cache_path = None

    file_names = []
    files = []

    for i in range(time_delta.days + 2):
        file_date = start_date + timedelta(days=i)
        file_name = str(file_date).replace('-', '') + '.nc'
        file_names.append(file_name)
        cache_file = cache_path + file_name
        ftp_file = path + file_name

        if os.path.isfile(cache_file):
            files.append(cache_file)
        else:
            print(f'Attempting to download {ftp_file}... ', end='')
            if cache and os.path.exists(cache_path):
                files.append(urlretrieve(ftp_file, cache_file)[0])
            else:
                files.append(urlretrieve(path + file_name)[0])
            print('done.')

    return files


def linear_interpolation_on_uniform_regular_grid(grid_info, point, *data):

    data_list = [None] * len(data)

    for dim in range(len(point)):

        x0, dx, xn = grid_info[dim]
        xi = point[dim]

        try:
            assert x0 <= xi <= xn, f'Point out of range in dim {dim} ({xi} is not in ({x0}, {xn})).'
        except TypeError as e:
            print(e, f'(in dim {dim}: x0 = {x0}, xi = {xi}, and xn = {xn})')

        index = (xi - x0) / dx
        index_floor = int(np.floor(index))
        index_diff = index - index_floor

        for i, data_ in enumerate(data):
            data_slice = data_[index_floor: index_floor + 2, ...]
            data_list[i] = (1 - index_diff) * data_slice[0, ...] + index_diff * data_slice[1, ...]

        data = tuple(data_list)

    if len(data) == 1:
        return data[0]

    else:
        return data


class Interpolate:

    def __init__(self, grid_vectors, *data, **kwargs):

        self.data = data
        self.grid_vectors = grid_vectors
        self.reference_time = kwargs.pop('reference_time', np.datetime64('1950-01-01T00:00'))
        self.time_units = kwargs.pop('time_units', 'h')
        self.grid_info = get_grid_info(grid_vectors, **kwargs)
        self.xarray_interp = kwargs.pop('xarray_interp', False)

    def interpolate(self, point):

        point_list = []

        for p in point:
            if isinstance(p, np.datetime64):
                p = (p - self.reference_time) / np.timedelta64(1, self.time_units)
            point_list.append(p)

        point = tuple(point_list)

        if self.xarray_interp:
            values = []
            for data in self.data:
                values.append(data.interp(time=point[0], latitude=point[1], longitude=point[2], assume_sorted=True))
            return tuple(values)

        else:
            return linear_interpolation_on_uniform_regular_grid(self.grid_info, point, *self.data)


def get_grid_info(grid_vectors, **kwargs):

    reference_time = kwargs.pop('reference_time', np.datetime64('1950-01-01T00:00'))
    time_units = kwargs.pop('time_units', 'h')

    grid_info = []

    for grid_vector in grid_vectors:

        if isinstance(grid_vector, xr.core.dataarray.DataArray):
            grid_vector = grid_vector.values

        if isinstance(grid_vector[0], np.datetime64):
            grid_vector = (grid_vector - reference_time) / np.timedelta64(1, time_units)

        grid_vector_min = np.min(grid_vector[0])
        grid_vector_max = np.max(grid_vector[-1])
        grid_vector_step = np.mean(np.diff(grid_vector))
        grid_info.append([grid_vector_min, grid_vector_step, grid_vector_max])

    return grid_info


def compute_magnitude(eastward_component, northward_component):

    return np.sqrt(eastward_component**2 + northward_component**2)


def compute_direction(eastward_component, northward_component):
    # computes angle clockwise from true north
    return np.rad2deg(np.arctan2(northward_component, eastward_component) - np.pi / 2)
