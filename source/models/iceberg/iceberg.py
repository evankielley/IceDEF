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
        
    Todo:
        * Add cross-sectional area info into get_shape_info function then remove get_cross_sectional areas

    """
    
    density = 900
    keel_shape = 'triangular'
    sail_shape = 'rectangular'
    
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
        self.length, self.width, self.sail_height = self.get_berg_dims(size)
        self.shape = shape
        self.shape_factor, self.height2draft_ratio, self.height, self.keel_depth = self.get_shape_info(shape, self.sail_height)
        self.mass = self.length*self.width*self.height*self.density
        self.bottom_area, self.top_area, self.keel_area, self.sail_area = self.get_cross_sectional_areas(self.length, self.width, self.sail_height, self.keel_depth)
        self.air_drag_coeff, self.water_drag_coeff, self.air_skin_drag_coeff, self.water_skin_drag_coeff = self.get_drag_coeffs() 
        
    
    def get_cross_sectional_areas(self, length, width, sail_height, keel_depth):
        bottom_area = 0
        top_area = length*width
        keel_area = keel_depth*length/2
        sail_area = sail_height*length
        return bottom_area, top_area, keel_area, sail_area
        
    
    def get_drag_coeffs(vary_air=True, vary_water=True):
        
        air_skin_drag_coeff = 2.5e-4
        water_skin_drag_coeff = 5.0e-4
        
        if vary_air:
            air_drag_coeff = np.random.uniform(0.5,2.5)
        else:
            air_drag_coeff = 1.5
            
        if vary_water:
            water_drag_coeff = np.random.uniform(0.5,2.5)
        else:
            water_drag_coeff = 1.5
            
            
        return air_drag_coeff, water_drag_coeff, air_skin_drag_coeff, water_skin_drag_coeff
    
    def get_shape_info(self, shape, sail_height):
        if shape == 'BLK' or 'TAB' or 'ISL' or 'RAD' or 'GEN':
            # no info for ISL, RAD, or GEN so assume BLK
            shape_factor = 0.5
            height2draft_ratio = 1/5
            keel_depth = sail_height/height2draft_ratio
            height = sail_height + keel_depth
        elif shape =='NTB':
            shape_factor = 0.41
            height2draft_ratio = 1/5
            keel_depth = sail_height/height2draft_ratio
            height = sail_height + keel_depth
        elif shape == 'DOM':
            shape_factor = 0.41
            height2draft_ratio = 1/4
            keel_depth = sail_height/height2draft_ratio
            height = sail_height + keel_depth
        elif shape == 'WDG':
            shape_factor = 0.33
            height2draft_ratio = 1/5
            keel_depth = sail_height/height2draft_ratio
            height = sail_height + keel_depth
        elif shape == 'PIN':
            shape_factor = 0.25
            height2draft_ratio = 1/2
            keel_depth = sail_height/height2draft_ratio
            height = sail_height + keel_depth
        elif shape == 'DD':
            shape_factor = 0.15
            height2draft_ratio = 1/1
            keel_depth = sail_height/height2draft_ratio
            height = sail_height + keel_depth
        else:
            print('Unknown shape {}'.format(shape))
            
        return shape_factor, height2draft_ratio, height, keel_depth
    
            
    def get_berg_dims(self, size, vary=True):
        # Size must be GR, BB, SM, MED, LG, VLG, or a list of [l, w, h]
        # See https://nsidc.org/data/g00807 for more info
        
        # Note: h is sail height
        
        if type(size) == list and len(size) == 3:
            l, w, h = size[0], size[1], size[2]
            
        elif type(size) == str and vary:
            if size == 'GR':
                l = np.random.uniform(0,5); w = np.random.uniform(0,5); h = np.random.uniform(0,1)
            elif size == 'BB':
                l = np.random.uniform(5,15); w = np.random.uniform(5,15); h = np.random.uniform(1,5)
            elif size == 'SM':
                l = np.random.uniform(15,60); w = np.random.uniform(15,60); h = np.random.uniform(5,15)
            elif size == 'MED':
                l = np.random.uniform(60,120); w = np.random.uniform(60,120); h = np.random.uniform(15,45)
            elif size == 'LG':
                l = np.random.uniform(120,200); w = np.random.uniform(120,200); h = np.random.uniform(45,75)
            elif size == 'VLG':
                l = np.random.uniform(200,500); w = np.random.uniform(200,500); h = np.random.uniform(75,150)
            elif size == 'GEN':
                # no info for GEN so I assume it's LG
                l = np.random.uniform(120,200); w = np.random.uniform(120,200); h = np.random.uniform(450,750)
            else:
                print('Unknown size class')
        
        elif type(size) == str and not vary:
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
                
        return l, w, h

def get_iip_df(season_year):
    # Choose iceberg year (2002 - 2015 available)
    # Note: Iceberg Season starts in November so many datasets include dates from year-1
    iip_url_base = 'ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G00807/' 
    iip_filename = 'IIP_{}IcebergSeason.csv'.format(season_year)
    iip_url = iip_url_base + iip_filename
    r = urllib.request.urlretrieve(iip_url)
    iip_df = pd.read_csv(r[0])
    return iip_df
    
def add_datetime_column(iip_df):
    iip_df['TIMESTAMP'] = pd.to_datetime(iip_df['SIGHTING_DATE'], format='%m/%d/%Y')
    iip_df = iip_df.loc[iip_df['SIGHTING_TIME'] >= 100]
    iip_df['TIMESTAMP'] += pd.to_timedelta(pd.to_datetime(iip_df['SIGHTING_TIME'], format='%H%M').dt.hour, unit='h')
    iip_df['TIMESTAMP'] += pd.to_timedelta(pd.to_datetime(iip_df['SIGHTING_TIME'], format='%H%M').dt.minute, unit='m')
    return iip_df

def get_time_dense_df(iip_df, max_hours):
    # max_hours is the desired max time between observations for individual icebergs (must be int)
    
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