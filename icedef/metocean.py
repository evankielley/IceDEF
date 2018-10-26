import os
import numpy as np
import xarray as xr
from urllib.request import urlretrieve
from datetime import date, timedelta


class Metocean:

    def __init__(self, date_bounds, **kwargs):
        self.date_bounds = date_bounds
        self.ocean_class = kwargs.pop('ocean_class', ECMWFOcean)
        self.atmosphere_class = kwargs.pop('atmosphere_class', NARRAtmosphere)
        self.ocean = self.ocean_class(date_bounds)
        self.atmosphere = self.atmosphere_class(date_bounds)

    def swap_ocean(self, ocean_class):
        self.ocean_class = ocean_class
        self.ocean = self.ocean_class(self.date_bounds)

    def swap_atmosphere(self, atmosphere_class):
        self.atmosphere_class = atmosphere_class
        self.atmosphere = self.atmosphere_class(self.date_bounds)


class ECMWFOcean:

    ID = "GLOBAL_ANALYSIS_FORECAST_PHY_001_024"
    PATH = 'ftp://data.munroelab.ca/pub/ECMWF/ocean/daily/'

    def __init__(self, date_bounds):

        self.dataset = xr.open_mfdataset(get_files(self.ID, self.PATH, date_bounds)).squeeze('depth').rename(
            {'uo': 'eastward_current_velocity', 'vo': 'northward_current_velocity'})

        self.eastward_current_velocities = xr.DataArray(data=self.dataset.eastward_current_velocity.values,
                                                        coords=[('time', self.dataset.time.values),
                                                                ('latitude', self.dataset.latitude.values),
                                                                ('longitude', self.dataset.longitude.values)],
                                                        attrs=self.dataset.eastward_current_velocity.attrs)
        self.northward_current_velocities = xr.DataArray(data=self.dataset.northward_current_velocity.values,
                                                         coords=[('time', self.dataset.time.values),
                                                                 ('latitude', self.dataset.latitude.values),
                                                                 ('longitude', self.dataset.longitude.values)],
                                                         attrs=self.dataset.eastward_current_velocity.attrs)


class ECMWFAtmosphere:

    ID = "WIND_GLO_WIND_L4_NRT_OBSERVATIONS_012_004"
    PATH = 'ftp://data.munroelab.ca/pub/ECMWF/atm/daily/'

    def __init__(self, date_bounds):

        self.dataset = xr.open_mfdataset(get_files(self.ID, self.PATH, date_bounds)).rename(
                                                        {'eastward_wind': 'eastward_wind_velocity',
                                                         'northward_wind': 'northward_wind_velocity'})

        self.eastward_wind_velocities = xr.DataArray(data=self.dataset.eastward_wind_velocity.values,
                                                     coords=[('time', self.dataset.time.values),
                                                             ('latitude', self.dataset.latitude.values),
                                                             ('longitude', self.dataset.longitude.values)],
                                                     attrs=self.dataset.eastward_wind_velocity.attrs)

        self.northward_wind_velocities = xr.DataArray(data=self.dataset.northward_wind_velocity.values,
                                                      coords=[('time', self.dataset.time.values),
                                                              ('latitude', self.dataset.latitude.values),
                                                              ('longitude', self.dataset.longitude.values)],
                                                      attrs=self.dataset.northward_wind_velocity.attrs)


class NARRAtmosphere:

    ID = "NCEP_North_American_Regional_Reanalysis_NARR"
    PATH = 'ftp://data.munroelab.ca/pub/NARR/atm/daily/'

    def __init__(self, date_bounds):

        self.dataset = xr.open_mfdataset(get_files(self.ID, self.PATH, date_bounds)).rename(
                                                        {'uwnd': 'eastward_wind_velocity',
                                                         'vwnd': 'northward_wind_velocity',
                                                         'lat': 'latitude',
                                                         'lon': 'longitude'})
        self.dataset['longitude'] = np.mod(self.dataset.longitude - 180, 360) - 180
        self.eastward_wind_velocities = xr.DataArray(data=self.dataset.eastward_wind_velocity.values,
                                                     coords=[('time', self.dataset.time.values),
                                                             ('latitude', self.dataset.latitude.values),
                                                             ('longitude', self.dataset.longitude.values)],
                                                     attrs=self.dataset.eastward_wind_velocity.attrs)

        self.northward_wind_velocities = xr.DataArray(data=self.dataset.northward_wind_velocity.values,
                                                      coords=[('time', self.dataset.time.values),
                                                              ('latitude', self.dataset.latitude.values),
                                                              ('longitude', self.dataset.longitude.values)],
                                                      attrs=self.dataset.northward_wind_velocity.attrs)


def get_files(id_, path, date_bounds, cache=True):

    start_date, end_date = date_bounds

    if isinstance(start_date, np.datetime64):
        start_date = start_date.astype(object)

    if isinstance(end_date, np.datetime64):
        end_date = end_date.astype(object)

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

    for i in range(time_delta.days + 1):
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

        assert x0 <= xi <= xn, f'Point out of range in dim {dim} ({xi} is not in ({x0}, {xn})).'

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

    def interpolate(self, point):

        point_list = []

        for p in point:
            if isinstance(p, np.datetime64):
                p = (p - self.reference_time) / np.timedelta64(1, self.time_units)
            point_list.append(p)

        point = tuple(point_list)

        return linear_interpolation_on_uniform_regular_grid(self.grid_info, point, *self.data)


def get_grid_info(grid_vectors, **kwargs):

    reference_time = kwargs.pop('reference_time', np.datetime64('1950-01-01T00:00'))
    time_units = kwargs.pop('time_units', 'h')

    grid_info = []

    for grid_vector in grid_vectors:

        if isinstance(grid_vector[0], np.datetime64):
            grid_vector = (grid_vector - reference_time) / np.timedelta64(1, time_units)

        if isinstance(grid_vector, xr.core.dataarray.DataArray):
            grid_vector = grid_vector.values

        grid_vector_min = grid_vector[0]
        grid_vector_max = grid_vector[-1]
        grid_vector_step = np.mean(np.diff(grid_vector))
        grid_info.append([grid_vector_min, grid_vector_step, grid_vector_max])

    return grid_info
