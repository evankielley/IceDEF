import urllib
import pandas as pd
import numpy as np


class Iceberg:
    
    density = 900
    keel_shape = 'triangular'
    sail_shape = 'rectangular'
    
    def __init__(self, id_num, datetimes, xvels, yvels, lats, lons, size, shape):
        self.id_num = id_num
        self.datetimes = datetimes
        self.xvels = xvels
        self.yvels = yvels
        self.lats = lats
        self.lons = lons
        self.size = size
        self.length, self.width, self.height = self.get_berg_dims(size)
        self.mass = self.length*self.width*self.height*self.density
        self.shape = shape
        self.shape_factor, self.height2draft_ratio, self.sail_height, self.keel_depth = self.get_shape_info(shape, self.height)
        self.bottom_area, self.top_area, self.sail_area, self.keel_area = self.get_cross_sectional_areas(self.length, self.width, self.sail_height, self.keel_depth)
        self.air_drag_coeff, self.water_drag_coeff, self.air_skin_drag_coeff, self.water_skin_drag_coeff = self.get_drag_coeffs()
        
    
    def get_cross_sectional_areas(self, length, width, sail_height, keel_depth):
        rhow = 1027.5
        bottom_area = 0
        top_area = length*width
        keel_area = keel_depth*length/2
        sail_area = sail_height*length
        return bottom_area, top_area, keel_area, sail_area
        
    
    def get_drag_coeffs(vary_air=True, vary_water=True):
        
        air_skin_drag_coeff = 2.5e-4
        water_skin_drag_coeff = 5.0e-4
        
        if vary_air:
            air_drag_coeff = np.random.uniform(0.01,2.5)
        else:
            air_drag_coeff = 1.2e-3
            
        if vary_water:
            water_drag_coeff = np.random.uniform(0.01,2.5)
        else:
            water_drag_coeff = 0.5
            
            
        return air_drag_coeff, water_drag_coeff, air_skin_drag_coeff, water_skin_drag_coeff
    
    def get_shape_info(self, shape, height):
        if shape == 'BLK' or 'TAB' or 'ISL' or 'RAD' or 'GEN':
            # no info for ISL, RAD, or GEN so assume BLK
            shape_factor = 0.5
            height2draft_ratio = 1/6
            sail_height = height2draft_ratio*height
            keel_depth = height - sail_height
        elif shape =='NTB':
            shape_factor = 0.41
            height2draft_ratio = 1/6
            sail_height = height2draft_ratio*height
            keel_depth = height - sail_height
        elif shape == 'DOM':
            shape_factor = 0.41
            height2draft_ratio = 1/5
            sail_height = height2draft_ratio*height
            keel_depth = height - sail_height
        elif shape == 'WDG':
            shape_factor = 0.33
            height2draft_ratio = 1/6
            sail_height = height2draft_ratio*height
            keel_depth = height - sail_height
        elif shape == 'PIN':
            shape_factor = 0.25
            height2draft_ratio = 1/3
            sail_height = height2draft_ratio*height
            keel_depth = height - sail_height
        elif shape == 'DD':
            shape_factor = 0.15
            height2draft_ratio = 1/2
            sail_height = height2draft_ratio*height
            keel_depth = height - sail_height
        else:
            print('Unknown shape {}'.format(shape))
        return shape_factor, height2draft_ratio, sail_height, keel_depth
    
            
    def get_berg_dims(self, size, vary=True):
        # Size must be GR, BB, SM, MED, LG, VLG, or a list of [l, w, h]
        # See https://nsidc.org/data/g00807 for more info
        
        if type(size) == list and len(size) == 3:
            l, w, h = size[0], size[1], size[2]
            
        elif type(size) == str and vary:
            if size == 'GR':
                l = np.random.uniform(0,5); w = np.random.uniform(0,5); h = np.random.uniform(0,10)
            elif size == 'BB':
                l = np.random.uniform(5,15); w = np.random.uniform(5,15); h = np.random.uniform(10,50)
            elif size == 'SM':
                l = np.random.uniform(15,60); w = np.random.uniform(15,60); h = np.random.uniform(50,150)
            elif size == 'MED':
                l = np.random.uniform(60,120); w = np.random.uniform(60,120); h = np.random.uniform(150,450)
            elif size == 'LG':
                l = np.random.uniform(120,200); w = np.random.uniform(120,200); h = np.random.uniform(450,750)
            elif size == 'VLG':
                l = np.random.uniform(200,500); w = np.random.uniform(200,500); h = np.random.uniform(750,1500)
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


def get_berg_df(season_year, chosen_track_ind):

    # Choose iceberg year (2002 - 2015 available)
    # Note: Iceberg Season starts in November so many datasets include dates from year-1
    iip_url_base = 'ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G00807/' 
    iip_filename = 'IIP_{}IcebergSeason.csv'.format(season_year)
    iip_url = iip_url_base + iip_filename
    r = urllib.request.urlretrieve(iip_url)
    iip_df = pd.read_csv(r[0], converters={'TIME':str})
    iip_df['TIMESTAMP'] = pd.to_datetime(iip_df['DATE'] + 'T' + iip_df['TIME'])


    # Choose the min number of observations for an eligible iceberg
    min_num_obs = 10
    eligible_bergs = np.asarray(
        iip_df['BERG_NUMBER'].value_counts()\
        .loc[iip_df['BERG_NUMBER'].value_counts() > min_num_obs].index)

    chosen_inds_arr = []

    for i in range(eligible_bergs.size):

        iip_berg_id = eligible_bergs[i]
        iip_berg_df = iip_df.loc[iip_df['BERG_NUMBER'] == iip_berg_id]
        
        ind0 = iip_berg_df.index.tolist()[0]
        indf = iip_berg_df.index.tolist()[-1]
        
        max_time_dif = np.timedelta64(24*60, 'm')
        
        chosen_inds = []

        for j in range(len(iip_berg_df)-1):

            time_dif = (iip_berg_df.TIMESTAMP.values[j+1] - \
                        iip_berg_df.TIMESTAMP.values[j]).astype('timedelta64[m]')
            
            if time_dif < max_time_dif:
                chosen_inds.append(j+ind0)

            elif len(chosen_inds) > 1:
                chosen_inds_arr.append(chosen_inds)
                chosen_inds = []
            else:
                chosen_inds = []

        if len(chosen_inds) > 1:
            chosen_inds_arr.append(chosen_inds)

    iip_berg_df = iip_df.loc[chosen_inds_arr[chosen_track_ind]].reset_index()

    return iip_berg_df
