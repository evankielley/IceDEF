"""This module can create an iceberg object and can extract a data frame from the IIP iceberg sighting database.
"""

import urllib
import pandas as pd
import numpy as np
import sys


class Iceberg():
    """Creates an iceberg object to be later used in drift simulation.

    Attributes:
        DIMS_BY_SIZE_CLASS (dict): iceberg dimension ranges according to IIP size classes
        RATIOS_BY_SHAPE_CLASS (dict): iceberg shape factors and height to draft ratios according to IIP shape classes
    """
    
    # dictionary values are of the form: [L_min, L_max, W_min, W_max, Hs_min, Hs_max]
    DIMS_BY_SIZE_CLASS = {'GR': [0, 5, 0, 5, 0, 1],
                          'BB': [5, 15, 5, 15, 1, 5],
                          'SM': [15, 60, 15, 60, 5, 15],
                          'MED': [60, 120, 60, 120, 15, 45],
                          'LG': [120, 200, 120, 200, 45, 75],
                          'VLG': [200, 400, 200, 400, 75, 150],
                          'GEN': [120, 200, 120, 200, 45, 75]}

    # dictionary values are of the form: [SF, H2D], where SF is shape factor and H2D is the height to draft ratio
    RATIOS_BY_SHAPE_CLASS = {'BLK': [0.5, 1/5], 
                             'TAB': [0.5, 1/5],
                             'ISL': [0.5, 1/5],
                             'RAD': [0.5, 1/5],
                             'GEN': [0.5, 1/5],
                             'NTB': [0.41, 1/5],
                             'DOM': [0.41, 1/4],
                             'WDG': [0.33, 1/5],
                             'PIN': [0.25, 1/3],
                             'DD': [0.15, 1/1]}
            

    
    def __init__(self, ID, t, x, y, vx, vy, size, shape):
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
            t (datetime.datetime): datetime of the iceberg
            x (float): iceberg longitude
            y (float): iceberg latitude
            vx (float): x-component of iceberg velocity (m/s)
            vy (float): y-component of iceberg velocity (m/s)
            size (str): size of the iceberg (can be GR, BB, MED, LG, VLG, or GEN)
            shape (str): shape of the iceberg. Can be BLK, TAB, ISL, GEN, RAD, NTB, DOM, WDG, PIN, or DD
        """
        
        self.keel_shape = 'triangular'
        self.sail_shape = 'rectangular'
        self.rho = 900  
        self.Cda = 1.5
        self.Cdw = 1.5
        self.Csda = 2.5e-4
        self.Csdw = 5.0e-4      
        
        self.ID = ID
        self.t = t
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.vx = vx
        self.ax = 0
        self.ay = 0
        self.size = size  # static -- do not change from init value, it has dep vars
        self.shape = shape  # static --do not change from init value, it has dep vars
        
        self.L, self.W, self.Hs = self.get_berg_dims()
        
        self.history = {'t': [], 'x': [], 'y': [], 'vx': [], 'vy': [], 'ax': [], 'ay': []}
        
        self.out_of_bounds = False
 
    
    @property
    def SF(self):
        return self.RATIOS_BY_SHAPE_CLASS[self.shape][0]
    
    @property
    def H2D(self):
        return self.RATIOS_BY_SHAPE_CLASS[self.shape][1]
    
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
            
            L_min, L_max, W_min, W_max, Hs_min, Hs_max = self.DIMS_BY_SIZE_CLASS[self.size]
                 
            L = (L_min + L_max)/2
            W = (W_min + W_max)/2
            Hs = (Hs_min + Hs_max)/2
        
        else:
            print("""Invalid size type. Either specify a list of [L, W, Hs] 
            or a str with an IIP size classification.""")
    
        return L, W, Hs

            
    def vary_L(self):
        L_min, L_max = self.DIMS_BY_SIZE_CLASS[self.size][0:2]
        self.L = np.random.uniform(L_min, L_max)
        
    def vary_W(self):
        W_min, W_max = self.DIMS_BY_SIZE_CLASS[self.size][2:4]
        self.W = np.random.uniform(W_min, W_max)

    def vary_Hs(self):
        Hs_min, Hs_max = self.DIMS_BY_SIZE_CLASS[self.size][4:6]
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
        
    def update_history(self):
        self.history['t'].append(self.t)
        self.history['x'].append(self.x)
        self.history['y'].append(self.y)
        self.history['vx'].append(self.vx)
        self.history['vy'].append(self.vy)
        self.history['ax'].append(self.ax)
        self.history['ay'].append(self.ay)

    def in_bounds(self, x_bounds, y_bounds):
        
        xmin = x_bounds[0]
        xmax = x_bounds[1]
        ymin = y_bounds[0]
        ymax = y_bounds[1]

        if not xmin < self.x < xmax:
            print('Iceberg out-of-bounds')
            return False

        elif not ymin < self.y < ymax:
            print('Iceberg out-of-bounds')
            return False

        else:
            return True
    
def clone_iceberg_state(berg):
    """This function clones the current state of an iceberg and returns the clone.
    
    Args:
        berg (icedef.iceberg.Iceberg): Iceberg object to be cloned.
        
    Returns:
        clone (icedef.iceberg.Iceberg): clone of the current state of the iceberg provided.
    """
    
    clone = Iceberg(berg.ID, berg.t, berg.x, berg.y, berg.vx, berg.vy, berg.size, berg.shape)
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


def get_dense_df(iip_season, max_hours, save=False):

    iip_df = add_datetime_column(get_iip_df(iip_season))

    iip_df = iip_df.sort_values(['ICEBERG_NUMBER', 'TIMESTAMP'], ascending=[True, True])

    iip_df = iip_df.reset_index(drop=True).drop(labels=['ICEBERG_YEAR', 'SIGHTING_DATE','SIGHTING_TIME','SIGHTING_METHOD'], axis=1)

    berg_nums = iip_df['ICEBERG_NUMBER'].tolist()
    berg_times = iip_df['TIMESTAMP'].tolist()

    berg_num = berg_nums[0]

    good_indices = []


    for i, row in iip_df.iterrows():

        if i+1 >= len(iip_df):
            break

        berg_num0 = berg_nums[i]
        berg_num1 = berg_nums[i+1]

        if berg_num0 == berg_num1:

            time0 = berg_times[i]
            time1 = berg_times[i+1]
            dtime = time1 - time0
            dt_hours = dtime.days*24 + dtime.seconds/3600

            if dt_hours < max_hours:
                good_indices.append(i)
                good_indices.append(i+1)

    good_indices = sorted(list(set(good_indices)))

    iip_df2 = iip_df[iip_df.index.isin(good_indices)]

    iip_df2['count'] = iip_df2.groupby('ICEBERG_NUMBER')['ICEBERG_NUMBER'].transform('count')

    iip_df2 = iip_df2.sort_values(['ICEBERG_NUMBER','TIMESTAMP'], ascending=[True, True])

    iip_df2 = iip_df2.reset_index(drop=True)

    berg_nums = iip_df2['ICEBERG_NUMBER'].tolist()
    berg_times = iip_df2['TIMESTAMP'].tolist()

    track_num = 0
    iip_df2['track_num'] = pd.Series(dtype=int)

    for i, row in iip_df2.iterrows():
        if i+2 > len(iip_df2):
            break
        berg_num0 = berg_nums[i]
        berg_num1 = berg_nums[i+1]
        time0 = berg_times[i]
        time1 = berg_times[i+1]
        dtime = time1 - time0
        dt_hours = dtime.days*24 + dtime.seconds/3600
        if berg_num0 == berg_num1 and dt_hours < max_hours:
            iip_df2.loc[i, 'track_num'] = track_num
            iip_df2.loc[i+1, 'track_num'] = track_num
        else:
            track_num += 1


    iip_df2 = iip_df2.sort_values(['track_num'], ascending=[True]).reset_index(drop=True)
    
    if save:
        iip_df2.to_csv(f'csvs/{iip_season}_max{max_hours}hr_tracks')

    return iip_df2