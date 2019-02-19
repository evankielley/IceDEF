import numpy as np
import xarray as xr
from icedef import statoil_arcticnet_data as sd


class TestCase:

    def __init__(self, beacon_id='204980', start_time=None, end_time=None):

        self.df = sd.get_beacon_df(sd.beacon_dir_path + '0' + beacon_id + '0' + '_2015.csv', start_time=start_time, end_time=end_time)
        self.ref_berg = sd.create_ref_berg_from_df(self.df)
        self.ref_times = self.ref_berg.history['time']
        self.ref_lats = xr.DataArray(self.ref_berg.history['latitude'], coords=[self.ref_times], dims=['time'])
        self.ref_lons = xr.DataArray(self.ref_berg.history['longitude'], coords=[self.ref_times], dims=['time'])
        self.start_time = np.datetime64(self.df['DataDate_UTC'].values[0])
        self.end_time = np.datetime64(self.df['DataDate_UTC'].values[len(self.df)-1])
        self.time_frame = self.start_time, self.end_time
        self.start_latitude = self.df.Latitude[0]
        self.start_longitude = self.df.Longitude[0]
        self.start_location = self.start_latitude, self.start_longitude
        self.start_velocity = sd.get_iceberg_velocity_from_dataframe(self.df, 0, 10)


class TestCaseA(TestCase):

    START_TIME = np.datetime64('2015-05-06T15:27:39')
    END_TIME = np.datetime64('2015-05-07T06:25:51')
    BEACON_ID = '90679'
    BEACON_CSV_FILE_PATH = sd.beacon_dir_path + '0' + BEACON_ID + '_2015.csv'

    def __init__(self):

        super().__init__(beacon_id=self.BEACON_ID, start_time=self.START_TIME, end_time=self.END_TIME)
