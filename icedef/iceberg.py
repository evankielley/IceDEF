"""This module can create an iceberg object and can extract a data frame from the IIP iceberg sighting database.
"""

import urllib
import pandas as pd
import numpy as np
import sys


class Iceberg():
    """Creates an iceberg object to be later used in drift simulation.

    Attributes:
        size_dim_dict (dict): iceberg dimension ranges according to IIP size classes
        shape_info_dict (dict): iceberg shape factors and height to draft ratios according to IIP shape classes
    """
    
    # dictionary values are of the form: [L_min, L_max, W_min, W_max, Hs_min, Hs_max]
    size_dim_dict = {'GR': [0, 5, 0, 5, 0, 1],
                     'BB': [5, 15, 5, 15, 1, 5],
                     'SM': [15, 60, 15, 60, 5, 15],
                     'MED': [60, 120, 60, 120, 15, 45],
                     'LG': [120, 200, 120, 200, 45, 75],
                     'VLG': [200, 400, 200, 400, 75, 150],
                     'GEN': [120, 200, 120, 200, 45, 75]}

    # dictionary values are of the form: [SF, H2D], where SF is shape factor and H2D is the height to draft ratio
    shape_info_dict = {'BLK': [0.5, 1/5], 
                       'TAB': [0.5, 1/5],
                       'ISL': [0.5, 1/5],
                       'RAD': [0.5, 1/5],
                       'GEN': [0.5, 1/5],
                       'NTB': [0.41, 1/5],
                       'DOM': [0.41, 1/4],
                       'WDG': [0.33, 1/5],
                       'PIN': [0.25, 1/3],
                       'DD': [0.15, 1/1]}
            

    
    def __init__(self, ID, T, X, Y, Vx, Vy, Ax, Ay, size, shape):
        """Instantiate iceberg object with necessary initial values.
        
        Attributes:
            keel_shape (str): shape of the iceberg keel
            sail_shape (str): shape of the iceberg sail
            rho (float): density of the iceberg
            Cda (float): air drag coefficient for the iceberg
            Cdw (float): water drag coefficient for the iceberg
            Csda (float): air skin drag coefficient for the iceberg
            Csdw (float): water skin drag coefficient for the iceberg

        Args:
            ID (int): iceberg ID number
            T (datetime.datetime): datetime of the iceberg
            Vx (float): x-component of iceberg velocity
            Vy (float): y-component of iceberg velocity
            Y (float): iceberg latitude
            X (float): iceberg longitude
            size (str): size of the iceberg (can be GR, BB, MED, LG, VLG, or GEN)
            shape (str): shape of the iceberg. Can be BLK, TAB, ISL, GEN, RAD, NTB, DOM, WDG, PIN, or DD
        """
        
        self.keel_shape = 'triangular'
        self.sail_shape = 'rectangular'
        self.rho = 900  
        self.Cda = 1.25
        self.Cdw = 1.25
        self.Csda = 2.5e-4
        self.Csdw = 5.0e-4      
        
        self.ID = ID
        self.T = T
        self.X = X
        self.Y = Y
        self.Vx = Vx
        self.Vy = Vy
        self.Ax = Ax
        self.Ay = Ay
        self.size = size  # static -- do not change from init value, it has dep vars
        self.shape = shape  # static --do not change from init value, it has dep vars
        
        self.L, self.W, self.Hs = self.get_berg_dims()
        
        self.history = {'T': [], 'X': [], 'Y': [], 'Vx': [], 'Vy': [], 'Ax': [], 'Ay': []}
 
    
    @property
    def SF(self):
        return self.shape_info_dict[self.shape][0]
    
    @property
    def H2D(self):
        return self.shape_info_dict[self.shape][1]
    
    @property
    def Hk(self):
        return self.Hs/self.H2D
    
    @property
    def H(self):
        return self.Hk + self.Hs
    
    @property
    def Ab(self):
        if self.keel_shape == 'rectangular':
            factor = 1
        elif self.keel_shape == 'triangular':
            factor = 0
        else:
            print('Invalid keel shape. Please choose triangular or rectangular.')
            
        return self.L*self.W*factor
    
    @property
    def Ak(self):
        if self.keel_shape == 'rectangular':
            factor = 1
        elif self.keel_shape == 'triangular':
            factor = 0.5
        else:
            print('Invalid keel shape. Please choose triangular or rectangular.')
            
        return 0.5*(self.L+self.W)*self.Hk*factor
    
    @property
    def At(self):
        if self.sail_shape == 'rectangular':
            factor = 1
        elif self.sail_shape == 'triangular':
            factor = 0
        else:
            print('Invalid sail shape. Please choose triangular or rectangular.')
            
        return self.L*self.W*factor
    
    @property
    def As(self):
        if self.sail_shape == 'rectangular':
            factor = 1
        elif self.sail_shape == 'triangular':
            factor = 0.5
        else:
            print('Invalid sail shape. Please choose triangular or rectangular.')
            
        # multiply by 0.5 to get average of L and W
        return 0.5*(self.L + self.W)*self.Hs*factor
    
    @property
    def M(self):
        #M = Mk + Ms; where Mk and Ms are the masses of the keel and sail, respectively.
        return self.rho*self.W*(self.Ak + self.As)
        
    
    def get_berg_dims(self):
        """This function returns numeric values for the dimensions of the iceberg (length, width, and sail height).
        
        Note:
            This function relies on there being a valid size attribute for the iceberg object.
            This can be a list of numeric values, [L, W, Hs] or a string representing a size classification code.
            See https://nsidc.org/data/g00807 for more information.
            Varying the dimensions can only happen if the iceberg is initialized with a size class (str) and not a list of dimensions.
            
        Returns:
            L (float): length of the iceberg (m)
            W (float): width of the iceberg (m)
            Hs (float): height of the sail of the iceberg (m)
        """
        
        
        if type(self.size) == list and len(self.size) == 3:
            L, W, Hs = self.size[0], self.size[1], self.size[2]
            
        elif type(self.size) == str:
            
            L_min, L_max, W_min, W_max, Hs_min, Hs_max = self.size_dim_dict[self.size]
                 
            L = (L_min + L_max)/2
            W = (W_min + W_max)/2
            Hs = (Hs_min + Hs_max)/2
        
        else:
            print("""Invalid size type. Either specify a list of [L, W, Hs] 
            or a str with an IIP size classification.""")
    
        return L, W, Hs

            
    def vary_L(self):
        L_min, L_max = self.size_dim_dict[self.size][0:2]
        self.L = np.random.uniform(L_min, L_max)
        
    def vary_W(self):
        W_min, W_max = self.size_dim_dict[self.size][2:4]
        self.W = np.random.uniform(W_min, W_max)

    def vary_Hs(self):
        Hs_min, Hs_max = self.size_dim_dict[self.size][4:6]
        self.Hs = np.random.uniform(Hs_min, Hs_max)
    
    def vary_all_dims(self):
        self.vary_L()
        self.vary_W()
        self.vary_Hs()
    
    def vary_Cda(self):
        self.Cda = np.random.uniform(0.5, 2.5)
        
    def vary_Cdw(self):
        self.Cdw = np.random.uniform(0.5, 2.5)
        
    def vary_all_drag_coeffs(self):
        self.vary_Cda()
        self.vary_Cdw()
       

    
def clone_iceberg_state(berg):
    """This function clones the current state of an iceberg and returns the clone.
    
    Args:
        berg (icedef.iceberg.Iceberg): Iceberg object to be cloned.
        
    Returns:
        clone (icedef.iceberg.Iceberg): clone of the current state of the iceberg provided.
    """
    
    clone = Iceberg(berg.ID, berg.T, berg.X, berg.Y, berg.Vx, berg.Vy, berg.Ax, berg.Ay, berg.size, berg.shape)
    return clone
    

def get_iip_df(season_year):
    """This function returns the IIP iceberg sighting database file for a particular year as a pandas data frame.
    
    Note:
       The iceberg seasons from 2002-2017 are available and the iceberg Season starts in November each year.
    
    Args:
        season_year (int): iceberg season year of the desired data frame
        
    Returns:
        iip_df (pandas.core.frame.DataFrame): IIP iceberg sighting database file for season_year specified
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
        iip_df (pandas.core.frame.DataFrame): IIP iceberg sighting dataframe
        
    Returns:
        iip_df (pandas.core.frame.DataFrame): IIP iceberg sighting dataframe with an added column, TIMESTAMP, with sighting T    
    """
    
    iip_df['TIMESTAMP'] = pd.to_datetime(iip_df['SIGHTING_DATE'], format='%m/%d/%Y')
    iip_df = iip_df.loc[iip_df['SIGHTING_TIME'] >= 100].copy()
    iip_df.loc[:, 'TIMESTAMP'] += pd.to_timedelta(pd.to_datetime(iip_df.loc[:,'SIGHTING_TIME'], format='%H%M').dt.hour, unit='h')
    iip_df.loc[:,'TIMESTAMP'] += pd.to_timedelta(pd.to_datetime(iip_df.loc[:,'SIGHTING_TIME'], format='%H%M').dt.minute, unit='m')
    
    return iip_df


def get_time_dense_df(iip_df, max_hours):
    """This function returns a new dataframe with only the rows of observations that are within the max time difference specified.
    
    Args:
        iip_df (pandas.core.frame.DataFrame): IIP iceberg sighting dataframe with an added column, TIMESTAMP, with sighting T
        max_hours (int): max number of hours desired between observations (hours)
        
    Returns:
        new_df (pandas.core.frame.DataFrame): IIP iceberg sighting dataframe with only rows of observations that are within the max hours specified
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


def get_iip_berg_df(iip_df, indices=range(3284, 3286)):
    """This function returns a dataframe from the desired rows of an IIP Iceberg Sighting dataframe.
    
    Args:
        iip_df (pandas.core.frame.DataFrame): dataframe from IIP Iceberg Sighting database
        indices (list or range): indices from the iip_df you wish to include in new dataframe
        
    Returns:
        iip_berg_df (pandas.core.frame.DataFrame): dataframe of just chosen rows from iip_df arg
    
    """
    
    if isinstance(indices, list):
        iip_berg_df = pd.DataFrame()
        for index in indices:
            iip_berg_df = iip_berg_df.append(iip_df.loc[iip_df.index == index])
            
    elif isinstance(indices, range):
        iip_berg_df = iip_df.loc[iip_df.index[indices]]
        
    else:
        print("indices argument type invalid -- indices must be either of type list or range")
        
    iip_berg_df = iip_berg_df.reset_index()

        
    return iip_berg_df


def get_iip_berg(iip_berg_df):
    """This function returns an Iceberg object (based off an IIP spreadsheet).
    
    Args:
        iip_berg_df (pandas.core.frame.DataFrame): dataframe of the IIP iceberg of interest
        
    Returns:
        iip_berg (icedef.iceberg.Iceberg): Iceberg object made from data in iip_berg_df arg 
    """
    
    
    ID = iip_berg_df['ICEBERG_NUMBER'].loc[0]
    T = iip_berg_df['TIMESTAMP'].dt.to_pydatetime()[0]
    X = iip_berg_df['SIGHTING_LONGITUDE'].loc[0]
    Y = iip_berg_df['SIGHTING_LATITUDE'].loc[0]
    Vx = 0
    Vy = 0
    Ax = 0
    Ay = 0
    size = iip_berg_df['SIZE'].loc[0]
    shape = iip_berg_df['SHAPE'].loc[0]
    
    iip_berg = Iceberg(ID, T, X, Y, Vx, Vy, Ax, Ay, size, shape)
    iip_berg.history['T'] = iip_berg_df['TIMESTAMP'].dt.to_pydatetime()
    iip_berg.history['X'] = iip_berg_df['SIGHTING_LONGITUDE'].loc[:].tolist()
    iip_berg.history['Y'] = iip_berg_df['SIGHTING_LATITUDE'].loc[:].tolist()
    
    return iip_berg