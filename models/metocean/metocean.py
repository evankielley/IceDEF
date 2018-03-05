import pandas as pd
import numpy as np
import netCDF4 as nc
import xarray as xr
import scipy.interpolate as interp
import datetime
from datetime import date, timedelta
import urllib


class Metocean(object):

    def __init__(self, x_min, x_max, y_min, y_max, t_min, t_max):
        self.x_min = x_min - abs(x_min-x_max) - self.xy_res
        self.x_max = x_max + abs(x_min-x_max) + self.xy_res
        self.y_min = y_min - abs(y_min-y_max) - self.xy_res
        self.y_max = y_max + abs(y_min-y_max) + self.xy_res
        
        if self.t_units == 'hours since 2000-01-01 00:00:00':
            self.t_min = (t_min - pd.Timestamp('2000-01-01')).days*24 + \
                         (t_min - pd.Timestamp('2000-01-01')).seconds/3600 - self.t_res
            self.t_max = (t_max - pd.Timestamp('2000-01-01')).days*24 + \
                         (t_max - pd.Timestamp('2000-01-01')).seconds/3600 + self.t_res
                
        elif self.t_units == 'numpy.datetime64':
            self.t_min = np.datetime64((t_min - pd.Timedelta('{}h'.format(self.t_res))))
            self.t_max = np.datetime64((t_max + pd.Timedelta('{}h'.format(self.t_res)))) 
            
        elif self.t_units == 'hours since 1950-01-01 00:00:00':
            self.t_min = (t_min - pd.Timestamp('1950-01-01')).days*24 + \
                         (t_min - pd.Timestamp('1950-01-01')).seconds/3600 - self.t_res
            self.t_max = (t_max - pd.Timestamp('1950-01-01')).days*24 + \
                         (t_max - pd.Timestamp('1950-01-01')).seconds/3600 + self.t_res
            self.t_min = nc.num2date(self.t_min, self.t_units, self.t_calendar)
            self.t_max = nc.num2date(self.t_max, self.t_units, self.t_calendar)
            
            self.year_min = str(self.t_min.year)
            self.year_max = str(self.t_max.year)
            self.month_min = str(self.t_min.month)
            self.month_max = str(self.t_max.month)
            self.day_min = str(self.t_min.day)
            self.day_max = str(self.t_max.day)

            if len(self.month_min) == 1:
                self.month_min = '0' + self.month_min
            if len(self.day_min) == 1:
                self.day_min = '0' + self.day_min
            if len(self.month_max) == 1:
                self.month_max = '0' + self.month_max
            if len(self.day_max) == 1:
                self.day_max = '0' + self.day_max
                
        elif self.t_units == 'hours since 1900-01-01 00:00:00.0 00:00':
            self.t_min = (t_min - pd.Timestamp('1900-01-01')).days*24 + \
                         (t_min - pd.Timestamp('1900-01-01')).seconds/3600 - self.t_res
            self.t_max = (t_max - pd.Timestamp('1900-01-01')).days*24 + \
                         (t_max - pd.Timestamp('1900-01-01')).seconds/3600 + self.t_res
            self.t_min = nc.num2date(self.t_min, self.t_units, self.t_calendar)
            self.t_max = nc.num2date(self.t_max, self.t_units, self.t_calendar)
            self.year_min = str(self.t_min.year)
            self.year_max = str(self.t_max.year)
            self.month_min = str(self.t_min.month)
            self.month_max = str(self.t_max.month)
            self.day_min = str(self.t_min.day)
            self.day_max = str(self.t_max.day)

            if len(self.month_min) == 1:
                self.month_min = '0' + self.month_min
            if len(self.day_min) == 1:
                self.day_min = '0' + self.day_min
            if len(self.month_max) == 1:
                self.month_max = '0' + self.month_max
            if len(self.day_max) == 1:
                self.day_max = '0' + self.day_max



class GLBv008(Metocean):
    url = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_56.3'
    xy_res = 0.08  # spatial resolution in degrees lat/lon
    t_res = 3  # temporal resolution in hours
    t_units = 'hours since 2000-01-01 00:00:00'
    def __init__(self, x_min, x_max, y_min, y_max, t_min, t_max):
        super().__init__(x_min, x_max, y_min, y_max, t_min, t_max)
        self.ds = xr.open_dataset(self.url, decode_times=False).sel(depth=0.0, 
                  lat = slice(self.y_min, self.y_max), 
                  lon = slice(self.x_min, self.x_max), 
                  time = slice(self.t_min, self.t_max))
        self.times = self.ds.variables['time'].values[:]
        self.datetimes = nc.num2date(testobj.ds.variables['time'].values[:],
                                     testobj.ds.time.units,
                                     testobj.ds.time.calendar)
        self.t2000 = self.times
        self.lats = np.asarray(self.ds.lat)
        self.lons = np.asarray(self.ds.lon)
        self.water_u = np.asarray(self.ds.water_u)
        self.water_v = np.asarray(self.ds.water_v)
        self.water_temp = np.asarray(self.ds.water_temp)
        self.water_u_interp = interp.RegularGridInterpolator((self.t2000, self.lats, self.lons), self.water_u)
        self.water_v_interp = interp.RegularGridInterpolator((self.t2000, self.lats, self.lons), self.water_v)
        self.water_temp_interp = interp.RegularGridInterpolator((self.t2000, self.lats, self.lons), self.water_temp)




class Navgem(Metocean):
    url = 'http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdNavgem05D10mWind_LonPM180'
    xy_res = 0.5  # spatial resolution in degrees lat/lon
    t_res = 6  # temporal resolution in hours
    t_units = 'numpy.datetime64'
    def __init__(self, x_min, x_max, y_min, y_max, t_min, t_max):
        super().__init__(x_min, x_max, y_min, y_max, t_min, t_max)
        self.ds = xr.open_dataset(self.url)
        self.ds = self.ds.sel(time = slice(self.t_min, self.t_max),
                  latitude = slice(self.y_min, self.y_max), 
                  longitude = slice(self.x_min, self.x_max))
        self.times = self.ds.variables['time']
        self.datetimes = []
        for i in (self.times - np.datetime64(0, 's'))/ np.timedelta64(1, 's'):
            self.datetimes.append(datetime.datetime.utcfromtimestamp(i))
        self.lats = np.asarray(self.ds.latitude)
        self.lons = np.asarray(self.ds.longitude)
        self.wind_u = np.asarray(self.ds.wnd_ucmp_height_above_ground)
        self.wind_v = np.asarray(self.ds.wnd_vcmp_height_above_ground)




class ECMWF_Ocean(Metocean):
    
    # product identifier: GLOBAL_ANALYSIS_FORECAST_PHY_001_024

    path = 'ftp://data.munroelab.ca/pub/ECMWF/ocean/daily/'
    xy_res = 1/12  # spatial resolution in degrees lat/lon
    t_res = 1  # temporal resolution in hours
    t_units = 'hours since 1950-01-01 00:00:00'
    t_calendar = 'standard'
    
    def __init__(self, x_min, x_max, y_min, y_max, t_min, t_max):
        
        super().__init__(x_min, x_max, y_min, y_max, t_min, t_max)
        
        d1 = date(self.t_min.year, self.t_min.month, self.t_min.day)  # start date
        d2 = date(self.t_max.year, self.t_max.month, self.t_max.day)  # end date
        delta = d2 - d1  # timedelta
            
        self.filenames = []
        self.files = []
        
        for i in range(delta.days + 1):
            new_date = d1 + timedelta(days=i)
            self.filenames.append(self.path + str(new_date).replace('-', '') + '.nc')
            self.files.append(urllib.request.urlretrieve(self.filenames[i])[0])
            
        if len(self.files) == 1:
            self.ds = nc.Dataset(self.files[0])
        elif len(self.files) > 1:
            self.ds = nc.MFDataset(self.files)
        else:
            print('Error in file list')

        self.times = nc.num2date(self.ds.variables['time'][:], self.t_units, self.t_calendar)
        self.datetimes = self.times
        self.t2000 = []
        for i in (self.datetimes - datetime.datetime(2000,1,1)):
            self.t2000.append(i.days*24 + i.seconds/3600)
        self.lats = self.ds.variables['latitude'][:]
        self.lons = self.ds.variables['longitude'][:]
        self.water_u = np.asarray(self.ds.variables['uo'][:,0,:,:])
        self.water_v = np.asarray(self.ds.variables['vo'][:,0,:,:])
        self.water_temp = np.asarray(self.ds.variables['thetao'][:,0,:,:])
        self.water_u_interp = interp.RegularGridInterpolator((self.t2000, self.lats, self.lons), self.water_u)
        self.water_v_interp = interp.RegularGridInterpolator((self.t2000, self.lats, self.lons), self.water_v)
        self.water_temp_interp = interp.RegularGridInterpolator((self.t2000, self.lats, self.lons), self.water_temp)




class ECMWF_Atm(Metocean):
    
    # product identifier: WIND_GLO_WIND_L4_NRT_OBSERVATIONS_012_004
    path = 'ftp://data.munroelab.ca/pub/ECMWF/atm/daily/'
    xy_res = 1/4  # spatial resolution in degrees lat/lon
    t_res = 6  # temporal resolution in hours
    t_units = 'hours since 1900-01-01 00:00:00.0 00:00'
    t_calendar = 'standard'
    
    def __init__(self, x_min, x_max, y_min, y_max, t_min, t_max):
        
        super().__init__(x_min, x_max, y_min, y_max, t_min, t_max)
        
        d1 = date(self.t_min.year, self.t_min.month, self.t_min.day)  # start date
        d2 = date(self.t_max.year, self.t_max.month, self.t_max.day)  # end date
        delta = d2 - d1  # timedelta
            
        self.filenames = []
        self.files = []
        
        for i in range(delta.days + 1):
            new_date = d1 + timedelta(days=i)
            self.filenames.append(self.path + 'sub' + str(new_date).replace('-', '') + '.nc')
            self.files.append(urllib.request.urlretrieve(self.filenames[i])[0])
            
        if len(self.files) == 1:
            self.ds = nc.Dataset(self.files[0])
        elif len(self.files) > 1:
            self.ds = nc.MFDataset(self.files)
        else:
            print('Error in file list')

        self.times = nc.num2date(self.ds.variables['time'][:], self.t_units, self.t_calendar)
        self.datetimes = self.times
        self.t2000 = []
        for i in (self.datetimes - datetime.datetime(2000,1,1)):
            self.t2000.append(i.days*24 + i.seconds/3600)
        self.lats = self.ds.variables['latitude'][:]
        self.lons = self.ds.variables['longitude'][:]
        self.wind_u = np.asarray(self.ds.variables['eastward_wind'][:,0,:,:])
        self.wind_v = np.asarray(self.ds.variables['northward_wind'][:,0,:,:])
        self.wind_u_interp = interp.RegularGridInterpolator((self.t2000, self.lats, self.lons), self.wind_u)
        self.wind_v_interp = interp.RegularGridInterpolator((self.t2000, self.lats, self.lons), self.wind_v)
