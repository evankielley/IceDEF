import numpy as np
import xarray as xr
from icedef import statoil_arcticnet_data as sd


class TestCase:

    def __init__(self, index_range=(1100, 1120)):
        self.df = sd.get_df(sd.dir_path + sd.csv_filenames[2])
        self.start_index, self.end_index = index_range
        self.ref_berg = sd.create_ref_berg_from_df(self.df, self.start_index, self.end_index)
        self.ref_times = self.ref_berg.history['time']
        self.ref_lats = xr.DataArray(self.ref_berg.history['latitude'], coords=[self.ref_times], dims=['time'])
        self.ref_lons = xr.DataArray(self.ref_berg.history['longitude'], coords=[self.ref_times], dims=['time'])
        self.start_time = np.datetime64(self.df.DataDate_UTC[self.start_index])
        self.end_time = np.datetime64(self.df.DataDate_UTC[self.end_index])
        self.start_latitude = self.df.Latitude[self.start_index]
        self.start_longitude = self.df.Longitude[self.start_index]
        self.start_velocity = sd.get_iceberg_velocity_from_dataframe(self.df, self.start_index, self.start_index + 10)
