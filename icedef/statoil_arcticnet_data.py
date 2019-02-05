"""This module works with the iceberg drift data collected by the Statoil-ArcticNet research cruise in 2015.
"""

from pandas import read_csv, to_datetime
from numpy import datetime64
from icedef import iceberg, tools

dir_path = 'http://icedef.munroelab.ca/data/StatoilArcticNet/drift_tracks/'
csv_filenames = ['0204980_2015.csv', '0505190_2015.csv', '0906790_2015.csv', '0907780_2015.csv']
kml_filenames = ['0204980_2015_ln.kml', '0505190_2015_ln.kml', '0906790_2015_ln.kml', '0907780_2015_ln.kml']
metadata_filename = 'MunroeMetadata.csv'
metadata_path = dir_path + metadata_filename
sample_data_path = dir_path + csv_filenames[0]


def get_df(path=sample_data_path):
    """This function returns a pandas dataframe of the beacon data at the specified path.

    :param path: location of the data file.
    :type path: str
    :return df: data file as a pandas dataframe.
    :rtype df: pandas.core.frame.DataFrame
    """
    
    df = read_csv(path)
    df.loc[:, 'DataDate_UTC'] = to_datetime(df['DataDate_UTC'])
    
    return df


def create_ref_berg_from_df(df, start_index, end_index):
    
    start_time = datetime64(df.DataDate_UTC[start_index])
    start_latitude = df.Latitude[start_index]
    start_longitude = df.Longitude[start_index]
    
    ref_berg = iceberg.quickstart(start_time, (start_latitude, start_longitude))

    for i in range(end_index - start_index + 2):

        if not df.DataDate_UTC[start_index + i] == df.DataDate_UTC[start_index + i + 1]:
            
            ref_berg.time = datetime64(df.DataDate_UTC[start_index + i])
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