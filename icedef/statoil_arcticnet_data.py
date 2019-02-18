"""This module works with the iceberg drift data collected by the Statoil-ArcticNet research cruise in 2015.
"""

from icedef import iceberg, tools
from urllib.request import urlretrieve
import xarray as xr
import numpy as np
import pandas as pd

root_dir_path = 'http://icedef.munroelab.ca/data/StatoilArcticNet/'
beacon_dir_path = root_dir_path + 'drift_tracks/'
beacon_csv_filenames = ['0204980_2015.csv', '0505190_2015.csv', '0906790_2015.csv', '0907780_2015.csv']
beacon_kml_filenames = ['0204980_2015_ln.kml', '0505190_2015_ln.kml', '0906790_2015_ln.kml', '0907780_2015_ln.kml']
beacon_metadata_filename = 'MunroeMetadata.csv'
beacon_metadata_path = beacon_dir_path + beacon_metadata_filename
sample_beacon_data_path = beacon_dir_path + beacon_csv_filenames[0]
avos_data_path = root_dir_path + 'AVOS_2015.csv'
adcp_data_path = root_dir_path + 'Leg1_1501_ADCP/an1501_os150bb.nc'


def get_df(path=sample_beacon_data_path):
    """This function returns a pandas dataframe of the beacon data at the specified path.

    :param path: location of the data file.
    :type path: str
    :return df: data file as a pandas dataframe.
    :rtype df: pandas.core.frame.DataFrame
    """
    
    df = pd.read_csv(path)
    df.loc[:, 'DataDate_UTC'] = pd.to_datetime(df['DataDate_UTC'])
    
    return df


def get_avos_df(path=avos_data_path, add_components=True):

    df = pd.read_csv(path, sep=' ; ')
    bad_rows = np.where((df['Longitude'] > 180) | (df['Longitude'] < -180) |
                        (df['Latitude'] > 180) | (df['Latitude'] < -180))[0]

    df = df.drop(df.index[bad_rows])
    df['Date'] = df['Date'].str.replace('/', '-')
    df['Date'] = pd.to_datetime(df['Date'])
    df['Wind speed'] = pd.to_numeric(df['Wind speed'])
    df['Wind dir'] = pd.to_numeric(df['Wind dir'])

    df.columns = ['time', 'latitude', 'longitude', 'direction', 'speed', 'air_temperature',
                  'water_temperature', 'eastward_velocity', 'northward_velocity']

    if add_components:

        df['eastward_velocity'] = -df['speed'] * np.sin(np.deg2rad(df['direction']))
        df['northward_velocity'] = -df['speed'] * np.cos(np.deg2rad(df['direction']))

    return df


def get_avos_ds(path=avos_data_path):

    df = get_avos_df(path)

    ds = df.set_index('time').to_xarray()

    return ds


def get_adcp_ds(path=adcp_data_path):

    file_loc, message = urlretrieve(path)

    ds = xr.open_dataset(file_loc).rename({
        'lon': 'longitude',
        'lat': 'latitude',
        'u': 'eastward_velocity',
        'v': 'northward_velocity',
        'uship': 'eastward_ship_velocity',
        'vship': 'northward_ship_velocity',
    })

    return ds


def create_ref_berg_from_df(df, start_index, end_index):
    
    start_time = np.datetime64(df.DataDate_UTC[start_index])
    start_latitude = df.Latitude[start_index]
    start_longitude = df.Longitude[start_index]
    
    ref_berg = iceberg.quickstart(start_time, (start_latitude, start_longitude))

    for i in range(end_index - start_index + 2):

        if not df.DataDate_UTC[start_index + i] == df.DataDate_UTC[start_index + i + 1]:
            
            ref_berg.time = np.datetime64(df.DataDate_UTC[start_index + i])
            ref_berg.latitude = df.Latitude[start_index + i]
            ref_berg.longitude = df.Longitude[start_index + i]
            ref_berg.update_history()
    
    return ref_berg


def get_iceberg_velocity_from_dataframe(df, start_index, end_index):
    
    dt = (df.DataDate_UTC[end_index] - df.DataDate_UTC[start_index]).total_seconds()  
    dlat = df.Latitude[end_index] - df.Latitude[start_index]
    dlon = df.Longitude[end_index] - df.Longitude[start_index]
    
    mid_lat = (df.Latitude[end_index] + df.Latitude[start_index]) / 2
    
    dy = tools.dlat_to_dy(dlat)
    dx = tools.dlon_to_dx(dlon, mid_lat)
    
    vx = dx/dt
    vy = dy/dt
    
    v = vx, vy
    
    return v