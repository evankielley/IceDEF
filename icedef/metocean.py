"""This module can instantiate objects which contain ECMWF ocean and atmospheric data for a particular time and space range.
"""

import scipy.interpolate as interp
import datetime
from datetime import date, timedelta
import urllib
import netCDF4 as nc
import numpy as np


class Metocean(object):
    """This class acts as a superclass that defines the spatial and temporal bounds for the data of its subclasses.
    """

    def __init__(self, x_min, x_max, y_min, y_max, t_min, t_max):
        """Instantiate a metocean object with spatial and temporal bounds.
        
        Args:
            x_min (float): minimum line of longitude for data region (degrees)
            x_max (float): maximum line of longitude for data region (degrees)
            y_min (float): minimum line of latitude for data region (degrees)
            y_max (float): maximum line of latitude for data region (degrees)
            t_min (datetime.datetime): minimum time for data time space
            t_max (datetime.datetime): maximum time for data time space
        """
        self.x_min = x_min - abs(x_min-x_max) - self.xy_res
        self.x_max = x_max + abs(x_min-x_max) + self.xy_res
        self.y_min = y_min - abs(y_min-y_max) - self.xy_res
        self.y_max = y_max + abs(y_min-y_max) + self.xy_res
        self.t_min = t_min - timedelta(hours = self.t_res)
        self.t_max = t_max + timedelta(hours = self.t_res)

    def convert_datetime2time(self, t, t_units, t_calendar, t_offset=0):
                                                
        dt += timedelta(hours = t_offset)
        t = nc.date2num(t, t_units, t_calendar)
        
        return t



class ECMWF_Ocean(Metocean):
    """This class creates an object which contains ocean data for surface current velocity and SST amongst other attributes.
    
    Note: 
        Product identifier: GLOBAL_ANALYSIS_FORECAST_PHY_001_024
    
    Args:
        Metocean (class): parent class
        
    Attributes:
        path (str): path to the directory of data files needed
        xy_res (float): spatial resolution of the ocean model (degrees)
        t_res (float): temporal resolution of the ocean model (hours)
        t_units (str): time units used in NetCDF data files
        t_calendar (str): time calendar used in NetCDF files
    """
    

    path = 'ftp://data.munroelab.ca/pub/ECMWF/ocean/daily/'
    xy_res = 1/12  # spatial resolution in degrees lat/lon
    t_res = 1  # temporal resolution in hours
    t_units = 'hours since 1950-01-01 00:00:00'
    t_calendar = 'standard'
    
    def __init__(self, x_min, x_max, y_min, y_max, t_min, t_max):
        
        super().__init__(x_min, x_max, y_min, y_max, t_min, t_max)
        
        self.filenames, self.files = self.get_filenames(self.t_min, self.t_max, self.path)
        self.ds = nc.MFDataset(self.files)
        self.times = self.ds.variables['time'][:]
        self.datetimes = nc.num2date(self.times, self.t_units, self.t_calendar)
        self.t1900 = nc.date2num(self.datetimes, 'hours since 1900-01-01 00:00:00.0 00:00', 'standard')
        self.t1950 = nc.date2num(self.datetimes, 'hours since 1950-01-01 00:00:00.0 00:00', 'standard')
        self.t2000 = nc.date2num(self.datetimes, 'hours since 2000-01-01 00:00:00.0 00:00', 'standard')
        self.lats = self.ds.variables['latitude'][:]
        self.lons = self.ds.variables['longitude'][:]
        self.UW = np.asarray(self.ds.variables['uo'][:,0,:,:])
        self.VW = np.asarray(self.ds.variables['vo'][:,0,:,:])
        self.SST = np.asarray(self.ds.variables['thetao'][:,0,:,:])
        self.iUW = interp.RegularGridInterpolator((self.times, self.lats, self.lons), self.UW)
        self.iVW = interp.RegularGridInterpolator((self.times, self.lats, self.lons), self.VW)
        self.iSST = interp.RegularGridInterpolator((self.times, self.lats, self.lons), self.SST)
        self.mean_u = np.mean(self.ds.variables['uo'][:,0,:,:].flatten())  # mean u
        self.mean_v = np.mean(self.ds.variables['vo'][:,0,:,:].flatten())  # mean v
        self.mean_u_mag = np.mean(abs(self.ds.variables['uo'][:,0,:,:].flatten()))  # mean magnitude of u
        self.mean_v_mag = np.mean(abs(self.ds.variables['vo'][:,0,:,:].flatten()))  # mean magnitude of v
        self.mean_mag = np.sqrt(self.mean_u_mag**2 + self.mean_v_mag**2)  # mean magnitude of resultant (m/s)
        self.mean_dir = np.arctan(self.mean_v/self.mean_u)  # mean direction (radians)
        self.std_u = np.std(self.ds.variables['uo'][:,0,:,:].flatten())  # std of u
        self.std_v = np.std(self.ds.variables['vo'][:,0,:,:].flatten())  # std of v
        self.std_u_mag = np.std(abs(self.ds.variables['uo'][:,0,:,:].flatten()))  # std of magnitude of u
        self.std_v_mag = np.std(abs(self.ds.variables['vo'][:,0,:,:].flatten()))  # std of magnitude of v
    
    def get_filenames(self, t_min, t_max, path):
        """This function returns the NetCDF files needed to access the desired ocean data (and their filenames)
        
        Args:
            t_min (datetime.datetime): minimum time for the data time space
            t_max (datetime.datetime): maximum time for the data time space
            path (str): path to the directory that contains the necessary data files
            
        Returns:
            filenames (list of str): list of the filenames of the files returned
            files (list): list of files accessed through server
        """

        d1 = date(self.t_min.year, self.t_min.month, self.t_min.day)  # start date
        d2 = date(self.t_max.year, self.t_max.month, self.t_max.day)  # end date
        delta = d2 - d1  # timedelta

        filenames = []
        files = []

        for i in range(delta.days + 1):
            new_date = d1 + timedelta(days=i)
            filenames.append(path + str(new_date).replace('-', '') + '.nc')
            files.append(urllib.request.urlretrieve(filenames[i])[0])

        return filenames, files



class ECMWF_Atm(Metocean):
    """This class creates an object which contains atmospheric data for 10 meter wind velocity amongst other attributes.
    
    Note: 
        Product identifier: WIND_GLO_WIND_L4_NRT_OBSERVATIONS_012_004
    
    Args:
        Metocean (class): parent class
        
    Attributes:
        path (str): path to the directory of data files needed
        xy_res (float): spatial resolution of the atmospheric model (degrees)
        t_res (float): temporal resolution of the atmospheric model (hours)
        t_units (str): time units used in NetCDF data files
        t_calendar (str): time calendar used in NetCDF files
    """
    
    path = 'ftp://data.munroelab.ca/pub/ECMWF/atm/daily/'
    xy_res = 1/4  # spatial resolution in degrees lat/lon
    t_res = 6  # temporal resolution in hours
    t_units = 'hours since 1900-01-01 00:00:00.0 00:00'
    t_calendar = 'standard'
    
    def __init__(self, x_min, x_max, y_min, y_max, t_min, t_max):
        
        super().__init__(x_min, x_max, y_min, y_max, t_min, t_max)
        
        self.filenames, self.files = self.get_filenames(self.t_min, self.t_max, self.path)
        self.ds = nc.MFDataset(self.files)
        self.times = self.ds.variables['time'][:]
        self.datetimes = nc.num2date(self.times, self.t_units, self.t_calendar)
        self.t1900 = nc.date2num(self.datetimes, 'hours since 1900-01-01 00:00:00.0 00:00', 'standard')
        self.t1950 = nc.date2num(self.datetimes, 'hours since 1950-01-01 00:00:00.0 00:00', 'standard')
        self.t2000 = nc.date2num(self.datetimes, 'hours since 2000-01-01 00:00:00.0 00:00', 'standard')
        self.lats = self.ds.variables['latitude'][:]
        self.lons = self.ds.variables['longitude'][:]
        self.UA = np.asarray(self.ds.variables['eastward_wind'][:,0,:,:])
        self.VA = np.asarray(self.ds.variables['northward_wind'][:,0,:,:])
        self.iUA = interp.RegularGridInterpolator((self.times, self.lats, self.lons), self.UA)
        self.iVA = interp.RegularGridInterpolator((self.times, self.lats, self.lons), self.VA)
        self.mean_u = np.mean(self.ds.variables['eastward_wind'][:,0,:,:].flatten())  # mean u
        self.mean_v = np.mean(self.ds.variables['northward_wind'][:,0,:,:].flatten())  # mean v
        self.mean_u_mag = np.mean(abs(self.ds.variables['eastward_wind'][:,0,:,:].flatten()))  # mean magnitude of u
        self.mean_v_mag = np.mean(abs(self.ds.variables['northward_wind'][:,0,:,:].flatten()))  # mean magnitude of v
        self.mean_mag = np.sqrt(self.mean_u_mag**2 + self.mean_v_mag**2)  # mean magnitude of resultant (m/s)
        self.mean_dir = np.arctan(self.mean_v/self.mean_u)  # mean direction (radians)
        self.std_u = np.std(self.ds.variables['eastward_wind'][:,0,:,:].flatten())  # std of u
        self.std_v = np.std(self.ds.variables['northward_wind'][:,0,:,:].flatten())  # std of v
        self.std_u_mag = np.std(abs(self.ds.variables['eastward_wind'][:,0,:,:].flatten()))  # std of magnitude of u
        self.std_v_mag = np.std(abs(self.ds.variables['northward_wind'][:,0,:,:].flatten()))  # std of magnitude of v
        
    def get_filenames(self, t_min, t_max, path):
        """This function returns the NetCDF files needed to access the desired ocean data (and their filenames)
        
        Args:
            t_min (datetime.datetime): minimum time for the data time space
            t_max (datetime.datetime): maximum time for the data time space
            path (str): path to the directory that contains the necessary data files
            
        Returns:
            filenames (list of str): list of the filenames of the files returned
            files (list): list of files accessed through server
        """

        d1 = date(self.t_min.year, self.t_min.month, self.t_min.day)  # start date
        d2 = date(self.t_max.year, self.t_max.month, self.t_max.day)  # end date
        delta = d2 - d1  # timedelta

        filenames = []
        files = []

        for i in range(delta.days + 1):
            new_date = d1 + timedelta(days=i)
            filenames.append(path + 'sub' + str(new_date).replace('-', '') + '.nc')
            files.append(urllib.request.urlretrieve(filenames[i])[0])

        return filenames, files