import os
import netCDF4 as nc
import xarray as xr
from urllib.request import urlretrieve
from datetime import date, timedelta


class ECMWFOcean:

    ID = "GLOBAL_ANALYSIS_FORECAST_PHY_001_024"
    PATH = 'ftp://data.munroelab.ca/pub/ECMWF/ocean/daily/'
    SPATIAL_RESOLUTION = 1/12  # spatial resolution in degrees lat/lon
    TIME_RESOLUTION = 1  # temporal resolution in hours
    TIME_UNITS = 'hours since 1950-01-01 00:00:00'
    TIME_CALENDAR = 'standard'

    def __init__(self, date_bounds):

        self.dataset = xr.open_mfdataset(get_files(self.ID, self.PATH,
                                                   date_bounds))
        self.eastward_current_velocities = xr.DataArray(
            data=self.dataset.uo.values[:, 0, :, :],
            coords=[('time', self.dataset.time.values),
                    ('latitude', self.dataset.latitude.values),
                    ('longitude', self.dataset.longitude.values)],
            attrs=self.dataset.uo.attrs)
        self.northward_current_velocities = xr.DataArray(
            data=self.dataset.vo.values[:, 0, :, :],
            coords=[('time', self.dataset.time.values),
                    ('latitude', self.dataset.latitude.values),
                    ('longitude', self.dataset.longitude.values)],
            attrs=self.dataset.vo.attrs)
        self.sea_surface_temperatures = xr.DataArray(
            data=self.dataset.thetao.values[:, 0, :, :],
            coords=[('time', self.dataset.time.values),
                    ('latitude', self.dataset.latitude.values),
                    ('longitude', self.dataset.longitude.values)],
            attrs=self.dataset.thetao.attrs)


class ECMWFAtmosphere:

    ID = "WIND_GLO_WIND_L4_NRT_OBSERVTIONS_012_004"
    PATH = 'ftp://data.munroelab.ca/pub/ECMWF/atm/daily/'
    SPATIAL_RESOLUTION = 1/4  # spatial resolution in degrees lat/lon
    TIME_RESOLUTION = 6  # temporal resolution in hours
    TIME_UNITS = 'hours since 1900-01-01 00:00:00.0 00:00'
    TIME_CALENDAR = 'standard'

    def __init__(self, date_bounds):

        self.dataset = xr.open_mfdataset(get_files(self.ID, self.PATH,
                                                   date_bounds))

        self.eastward_wind_velocities = xr.DataArray(
            data=self.dataset.eastward_wind.values,
            coords=[('time', self.dataset.time.values),
                    ('latitude', self.dataset.latitude.values),
                    ('longitude', self.dataset.longitude.values)],
            attrs=self.dataset.eastward_wind.attrs)

        self.northward_wind_velocities = xr.DataArray(
            data=self.dataset.northward_wind.values[:, 0, :, :],
            coords=[('time', self.dataset.time.values),
                    ('latitude', self.dataset.latitude.values),
                    ('longitude', self.dataset.longitude.values)],
            attrs=self.dataset.northward_wind.attrs)


class NARRAtmosphere:

    ID = "NCEP_North_American_Regional_Reanalysis_NARR"
    PATH = 'ftp://data.munroelab.ca/pub/NARR/atm/daily/'
    SPATIAL_RESOLUTION = 1/4  # spatial resolution in degrees lat/lon
    TIME_RESOLUTION = 3  # temporal resolution in hours
    TIME_UNITS = 'hours since 1800-1-1 00:00:0.0'
    TIME_CALENDAR = 'standard'

    def __init__(self, date_bounds):

        self.dataset = xr.open_mfdataset(get_files(self.ID, self.PATH,
                                                   date_bounds))
        self.eastward_wind_velocities = xr.DataArray(
            data=self.dataset.uwnd.values,
            coords=[('time', self.dataset.time.values),
                    ('latitude', self.dataset.lat.values),
                    ('longitude', self.dataset.lon.values)],
            attrs=self.dataset.uwnd.attrs)
        self.northward_wind_velocities = xr.DataArray(
            data=self.dataset.uwnd.values,
            coords=[('time', self.dataset.time.values),
                    ('latitude', self.dataset.lat.values),
                    ('longitude', self.dataset.lon.values)],
            attrs=self.dataset.uwnd.attrs)


def get_files(id_, path, date_bounds, cache=True):

    start_date, end_date = date_bounds
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

    filenames = []
    files = []

    for i in range(time_delta.days + 1):
        file_date = start_date + timedelta(days=i)
        filename = str(file_date).replace('-', '') + '.nc'
        filenames.append(filename)
        cache_file = cache_path + filename

        if os.path.isfile(cache_file):
            files.append(cache_file)
        else:
            if cache and os.path.exists(cache_path):
                files.append(urlretrieve(path +
                                         filename, cache_path + filename)[
                                 0])
            else:
                files.append(urlretrieve(path + filename)[0])

    return files


def interpolate(data_array, point, method='linear'):
    interpolated_value = data_array.interp(coords={'time': point[0],
                                            'latitude': point[1],
                                            'longitude': point[2]},
                                    assume_sorted=True,
                                    method=method)
    return interpolated_value
