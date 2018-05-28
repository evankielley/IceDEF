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
    
    def __init__(self, ID, T, X, Y, Vx, Vy, size, shape):
        """Instantiate iceberg object with necessary initial values.

        Note:
            All parameters that are lists must align.

        Args:
            ID (int): iceberg ID number
            T (datetime): datetime of the iceberg
            Vx (float): x-component of iceberg velocity
            Vy (float): y-component of iceberg velocity
            Y (float): iceberg latitude
            X (float): iceberg longitude
            size (str): size of the iceberg (can be GR, BB, MED, LG, VLG, or GEN)
            shape (str): shape of the iceberg. Can be BLK, TAB, ISL, GEN, RAD, NTB, DOM, WDG, PIN, or DD
        """
        self.ID = ID
        self.T = T
        self.Vx = Vx
        self.Vy = Vy
        self.Y = Y
        self.X = X
        self.size = size
        self.shape = shape
        self.length, self.width, self.sail_height = self.get_berg_dims()
        self.shape_factor, self.height2draft_ratio, self.height, self.keel_depth, self.bottom_area, self.top_area, self.keel_area, self.sail_area, self.mass = self.get_shape_info()
        self.history = {'T': [], 'X': [], 'Y': [], 'Vx': [], 'Vy': []}
    
    def get_shape_info(self):
        """This function returns information pertaining to the shape and dimensions of the iceberg.
        
        Note:
            All units are SI. Dimensions are in meters, areas are meters squared, and masses are kilograms.
        
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
        
        mass = self.length*self.width*height*Iceberg.density
            
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
    
    
def clone_iceberg_state(berg):
    """This function clones the current state of an iceberg and returns the clone.
    
    Args:
        berg (icedef.iceberg.Iceberg): Iceberg object to be cloned.
        
    Returns:
        clone (icedef.iceberg.Iceberg): clone of the current state of the iceberg provided.
    """
    
    clone = Iceberg(berg.ID, berg.T, berg.X, berg.Y, berg.Vx, berg.Vy, berg.size, berg.shape)
    return clone
    

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
    This function will take information from the SIGHTING_DATE and SIGHTING_TIME columns and convert into T in a new column.
    
    Args:
        iip_df (Dataframe): IIP iceberg sighting dataframe
        
    Returns:
        iip_df (Dataframe): IIP iceberg sighting dataframe with an added column, TIMESTAMP, with sighting T    
    """
    
    iip_df['TIMESTAMP'] = pd.to_datetime(iip_df['SIGHTING_DATE'], format='%m/%d/%Y')
    iip_df = iip_df.loc[iip_df['SIGHTING_TIME'] >= 100]
    iip_df['TIMESTAMP'] += pd.to_timedelta(pd.to_datetime(iip_df['SIGHTING_TIME'], format='%H%M').dt.hour, unit='h')
    iip_df['TIMESTAMP'] += pd.to_timedelta(pd.to_datetime(iip_df['SIGHTING_TIME'], format='%H%M').dt.minute, unit='m')
    return iip_df


def get_time_dense_df(iip_df, max_hours):
    """This function returns a new dataframe with only the rows of observations that are within the max time difference specified.
    
    Args:
        iip_df (Dataframe): IIP iceberg sighting dataframe with an added column, TIMESTAMP, with sighting T
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


def get_iip_iceberg(iip_season=2015, method='index', identifier=range(3283, 3285)):
    """This function returns an Iceberg object (based off an IIP spreadsheet) and its associated DataFrame.
    
    Args:
        iip_season (int): iceberg season (year) of the IIP Iceberg Sighting data of interest
        identifier (range or int): can be either IIP spreadsheet indices as a range or an Iceberg ID as an int.
        method (str): method to be used to retrieve iceberg data from spreadsheet, can be either 'index' or 'ID'
        
    Returns:
        iip_berg_df (pandas.core.frame.DataFrame): DataFrame of the iceberg specified.
        iip_berg (icedef.iceberg.Iceberg): Iceberg object made from IIP spreadsheet data specified. 
    """
    
    iip_df = get_iip_df(iip_season)
    iip_df = add_datetime_column(iip_df)
    
    if method == 'index':
        if isinstance(identifier, range):
            iip_berg_df = iip_df.loc[iip_df.index[identifier]].reset_index()
        else:
            print('Invalid identifier for index method. Identifier should be range')
            
    elif method == 'ID':
        if isinstance(identifier, int):
            iip_berg_df = iip_df.loc[iip_df['ICEBERG_NUMBER'] == identifier].reset_index()
        else:
            print('Invalid identifier for ID method. Identifier should be int')
            
    ID = iip_berg_df['ICEBERG_NUMBER'].loc[0]
    T = iip_berg_df['TIMESTAMP'].dt.to_pydatetime()[0]
    X = iip_berg_df['SIGHTING_LONGITUDE'].loc[0]
    Y = iip_berg_df['SIGHTING_LATITUDE'].loc[0]
    Vx = 0
    Vy = 0
    size = iip_berg_df['SIZE'].loc[0]
    shape = iip_berg_df['SHAPE'].loc[0]
    
    iip_berg = Iceberg(ID, T, X, Y, Vx, Vy, size, shape)
    iip_berg.history['T'] = iip_berg_df['TIMESTAMP'].dt.to_pydatetime()
    iip_berg.history['X'] = iip_berg_df['SIGHTING_LONGITUDE'].loc[:].tolist()
    iip_berg.history['Y'] = iip_berg_df['SIGHTING_LATITUDE'].loc[:].tolist()
    
    return iip_berg_df, iip_berg
