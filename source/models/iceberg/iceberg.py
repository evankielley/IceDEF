"""This module can create an iceberg object and can extract a data frame from the IIP iceberg sighting database.
"""

import urllib
import pandas as pd
import numpy as np


class Iceberg:
    """Creates an iceberg object to be later used in drift simulation.

    Attributes:
        density (float): density of the iceberg
        keel_shape (str): shape of the iceberg keel
        sail_shape (str): shape of the iceberg sail
        air_drag_coeff (float): air drag coefficient for the iceberg
        water_drag_coeff (float): water drag coefficient for the iceberg
        air_skin_drag_coeff (float): air skin drag coefficient for the iceberg
        water_skin_drag_coeff (float): water skin drag coefficient for the iceberg
    """
    
    density = 900
    keel_shape = 'triangular'
    sail_shape = 'rectangular'
    air_drag_coeff = 1.25
    water_drag_coeff = 1.25
    air_skin_drag_coeff = 2.5e-4
    water_skin_drag_coeff = 5.0e-4
    
    def __init__(self, id_num, datetimes, xvels, yvels, lats, lons, size, shape):
        """Instantiate iceberg object with necessary initial values.

        Note:
            All parameters that are lists must align.

        Args:
            id_num (int): iceberg ID number
            datetimes (list of datetime.datetime): list of iceberg's datetimes
            xvels (list of float): list of x-components of iceberg velocites
            yvels (list of float): list of y-components of iceberg velocites
            lats (list of float): list of iceberg latitudes
            lons (list of float): list of iceberg longitudes
            size (str): size of the iceberg (can be GR, BB, MED, LG, VLG, or GEN)
            shape (str): shape of the iceberg. Can be BLK, TAB, ISL, GEN, RAD, NTB, DOM, WDG, PIN, or DD
        """
        
        self.id_num = id_num
        self.datetimes = datetimes
        self.xvels = xvels
        self.yvels = yvels
        self.lats = lats
        self.lons = lons
        self.size = size
        self.shape = shape
        self.length, self.width, self.sail_height = self.get_berg_dims()
        self.shape_factor, self.height2draft_ratio, self.height, self.keel_depth, self.bottom_area, self.top_area, self.keel_area, self.sail_area, self.mass = self.get_shape_info()
       
    
    def get_shape_info(self):
        """This function returns information pertaining to the shape and dimensions of the iceberg.
        
        Note:
            All units are ANSI. Dimensions are in meters, areas are meters squared, and masses are kilograms.
        
        Returns:
            shape_factor (float): constant that is specified according to iceberg shape classification
            height2draft_ratio (float): constant that dictates the ratio of the iceberg above water versus below
            height (float): total height of the iceberg (sail_height + keel_depth) (meters)
            keel_depth (float): depth of the iceberg keel (meters)
            bottom_area (float): area of the bottom face of the iceberg (meters squared)
            top_area (float): area of the top face of the iceberg (meters squared)
            keel_area (float): area of the keel of the iceberg (meters squared)
            sail_area (float): area of the sail of the iceberg (meters squared)
            mass (float): mass of the iceberg (kilograms)
        """
        
        shape = self.shape
        
        if shape == 'BLK' or 'TAB' or 'ISL' or 'RAD' or 'GEN':
            # no info for ISL, RAD, or GEN so assume BLK
            shape_factor = 0.5
            height2draft_ratio = 1/5
            keel_depth = self.sail_height/height2draft_ratio
            height = self.sail_height + keel_depth
        elif shape =='NTB':
            shape_factor = 0.41
            height2draft_ratio = 1/5
            keel_depth = self.sail_height/height2draft_ratio
            height = self.sail_height + keel_depth
        elif shape == 'DOM':
            shape_factor = 0.41
            height2draft_ratio = 1/4
            keel_depth = self.sail_height/height2draft_ratio
            height = self.sail_height + keel_depth
        elif shape == 'WDG':
            shape_factor = 0.33
            height2draft_ratio = 1/5
            keel_depth = self.sail_height/height2draft_ratio
            height = self.sail_height + keel_depth
        elif shape == 'PIN':
            shape_factor = 0.25
            height2draft_ratio = 1/2
            keel_depth = self.sail_height/height2draft_ratio
            height = self.sail_height + keel_depth
        elif shape == 'DD':
            shape_factor = 0.15
            height2draft_ratio = 1/1
            keel_depth = self.sail_height/height2draft_ratio
            height = self.sail_height + keel_depth
        else:
            print('Unknown shape {}'.format(shape))
            
        bottom_area = 0
        top_area = self.length*self.width
        keel_area = keel_depth*self.length/2
        sail_area = self.sail_height*self.length
        
        mass = self.length*self.width*self.height*Iceberg.density
            
        return shape_factor, height2draft_ratio, height, keel_depth, bottom_area, top_area, keel_area, sail_area, mass
    
            
    def get_berg_dims(self):
        """This function returns numeric values for the dimensions of the iceberg (length, width, and sail height).
        
        Note:
            This function relies on there being a valid size attribute for the iceberg object.
            This can be a list of numeric values, [l, w, sail_height] or a string representing a size classification code.
            See https://nsidc.org/data/g00807 for more information.
            
        Returns:
            l (float): length of the iceberg (meters)
            w (float): width of the iceberg (meters)
            h (float): height of the sail of the iceberg (meters)
        """
        
        size = self.size
        
        if type(size) == list and len(size) == 3:
            l, w, h = size[0], size[1], size[2]
            
        elif type(size) == str:
            if size == 'GR':
                l = (0+5)/2; w = (0+5)/2; h = (0+1)/2*10
            elif size == 'BB':
                l = (5+15)/2; w = (5+15)/2; h = (1+5)/2*10        
            elif size == 'SM':
                l = (15+60)/2; w = (15+60)/2; h = (5+15)/2*10        
            elif size == 'MED':
                l = (60+120)/2; w = (60+120)/2; h = (15+45)/2*10               
            elif size == 'LG':
                l = (120+200)/2; w = (120+200)/2; h = (45+75)/2*10                
            elif size == 'VLG':
                # Sizes have no listed upper bound
                l = (200+500/2)/2; w = (200+500/2)/2; h = (75+75/2)/2*10     
            # This info for GEN is wrong!
            elif size == 'GEN':
                l = (120+200)/2; w = (120+200)/2; h = (45+75)/2*10            
            else:
                print('unknown size class')
                l = None; w = None; h = None               
        else:
            l, w, h = None, None, None
            print('Invalid size entry. Please enter a list of floats [l, w, h] or a valid size class as a string')
                
        return l, w, h
    

def get_iip_df(season_year):
    """This function returns the IIP iceberg sighting database file for a particular year as a pandas data frame.
    
    Note:
       The iceberg seasons from 2002-2017 are available and the iceberg Season starts in November each year.
    
    Args:
        season_year (int): iceberg season year of the desired data frame
        
    Returns:
        iip_df (Dataframe): IIP iceberg sighting database file for season_year specified
    """
    
    iip_url_base = 'ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G00807/' 
    iip_filename = 'IIP_{}IcebergSeason.csv'.format(season_year)
    iip_url = iip_url_base + iip_filename
    r = urllib.request.urlretrieve(iip_url)
    iip_df = pd.read_csv(r[0])
    return iip_df
    
    
def add_datetime_column(iip_df):
    """This function adds a column of datetimes to an IIP iceberg sighting database data frame.
    
    The arg iip_df should be obtained through the get_iip_df function. 
    This function will take information from the SIGHTING_DATE and SIGHTING_TIME columns and convert into datetimes in a new column.
    
    Args:
        iip_df (Dataframe): IIP iceberg sighting dataframe
        
    Returns:
        iip_df (Dataframe): IIP iceberg sighting dataframe with an added column, TIMESTAMP, with sighting datetimes    
    """
    
    iip_df['TIMESTAMP'] = pd.to_datetime(iip_df['SIGHTING_DATE'], format='%m/%d/%Y')
    iip_df = iip_df.loc[iip_df['SIGHTING_TIME'] >= 100]
    iip_df['TIMESTAMP'] += pd.to_timedelta(pd.to_datetime(iip_df['SIGHTING_TIME'], format='%H%M').dt.hour, unit='h')
    iip_df['TIMESTAMP'] += pd.to_timedelta(pd.to_datetime(iip_df['SIGHTING_TIME'], format='%H%M').dt.minute, unit='m')
    return iip_df


def get_time_dense_df(iip_df, max_hours):
    """This function returns a new dataframe with only the rows of observations that are within the max time difference specified.
    
    Args:
        iip_df (Dataframe): IIP iceberg sighting dataframe with an added column, TIMESTAMP, with sighting datetimes
        max_hours (int): max number of hours desired between observations (hours)
        
    Returns:
        new_df (Dataframe): IIP iceberg sighting dataframe with only rows of observations that are within the max hours specified
    """
    
    max_timedelta = np.timedelta64(60*max_hours, 'm')
    new_df = pd.DataFrame()

    for iceberg_number in iip_df['ICEBERG_NUMBER'].unique():
        tmp_df = iip_df.loc[iip_df['ICEBERG_NUMBER'] == iceberg_number].reset_index(drop=True)
        for i in range(len(tmp_df) - 1):
            if (tmp_df['TIMESTAMP'][i+1] - tmp_df['TIMESTAMP'][i]) < max_timedelta:
                new_df = new_df.append(tmp_df.loc[i])
                new_df = new_df.append(tmp_df.loc[i+1])
    
    new_df['ICEBERG_NUMBER'] = new_df['ICEBERG_NUMBER'].astype(int)
    new_df['ICEBERG_YEAR'] = new_df['ICEBERG_YEAR'].astype(int)
    new_df['SIGHTING_TIME'] = new_df['SIGHTING_TIME'].astype(int)

    new_df = new_df.drop_duplicates().reset_index(drop=True)
                
    return new_df
