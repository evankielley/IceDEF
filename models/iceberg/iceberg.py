import urllib
import pandas as pd
import numpy as np


class Iceberg:
    
    def __init__(self, id_num, times, t2000, lats, lons, size):
        self.id_num = id_num
        self.times = times
        self.t2000 = t2000
        self.times_ref0 = self.get_times_ref0(self.times)
        self.lats = lats
        self.lons = lons
        
        if type(size) == str:
            self.length, self.width, self.height = self.get_berg_dims(size)
        elif type(size) == list and len(size) == 3:
            self.length, self.width, self.height = size[0], size[1], size[2]
        else:
            print('Invalid size argument')
            
    def get_times_ref0(self, times):
        times_ref0 = []
        t_offset = times[0].hour + times[0].minute/60  # from start of day
        for i in range(len(times)):
            times_ref0.append(round((times[i] - times[0]).days*24 + \
                                    (times[i] - times[0]).seconds/3600 + t_offset, 1))
        return times_ref0
    
    def set_times_ref0(self, times):
        self.times_ref0 = self.get_times_ref0(times)
            
    def get_berg_dims(self, size):
        # Size must be GR, BB, SM, MED, LG, VLG
        # See https://nsidc.org/data/g00807 for more info
        if size == 'GR':
            l = (0+5)/2; w = (0+5)/2; h = (0+1)/2*10
        elif size == 'BB':
            l = (5+15)/2; w = (5+15)/2; h = (1+5)/2*10        
        elif size == 'SM':
            l = (15+60)/2; w = (15+60)/2; h = (5+15)/2*10        
        elif size == 'MED':
            l = (60+120)/2; w = (60+120)/2; h = (15+45)/2*10               
        elif size == 'LG':
            l = (120)/2; w = (120)/2; h = (45+75)/2*10                
        elif size == 'VLG':
            # Sizes have no listed upper bound
            l = (200+200/2)/2; w = (200+200/2)/2; h = (75+75/2)/2*10     
        # This info for GEN is wrong!
        elif size == 'GEN':
            l = (120)/2; w = (120)/2; h = (45+75)/2*10            
        else:
            print('unknown size class')
            l = None; w = None; h = None
        return l, w, h


def get_berg_df(chosen_track_ind):

    # Choose iceberg year (2002 - 2015 available)
    # Note: Iceberg Season starts in November so many datasets include dates from year-1
    season_year = 2015
    iip_url_base = 'ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G00807/' 
    iip_filename = 'IIP_{}IcebergSeason.csv'.format(season_year)
    iip_url = iip_url_base + iip_filename
    r = urllib.request.urlretrieve(iip_url)
    iip_df = pd.read_csv(r[0], converters={'TIME':str})
    iip_df['DATETIME'] = pd.to_datetime(iip_df['DATE'] + 'T' + iip_df['TIME'])
    iip_df['t2000'] = round((iip_df.DATETIME - pd.to_datetime('2000-01-01')).dt.days*24 + \
                      (iip_df.DATETIME - pd.to_datetime('2000-01-01')).dt.seconds/3600, 2)

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

            time_dif = (iip_berg_df.DATETIME.values[j+1] - \
                        iip_berg_df.DATETIME.values[j]).astype('timedelta64[m]')
            
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
    iip_berg_df['t2000'] = round((iip_berg_df.DATETIME - pd.to_datetime('2000-01-01')).dt.days*24 + \
                                 (iip_berg_df.DATETIME - pd.to_datetime('2000-01-01')).dt.seconds/3600, 2)

    return iip_berg_df
