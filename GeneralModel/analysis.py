"""
This module does analysis on the output of an Iceberg simulation.
"""

import matplotlib.pyplot as plt                                         
import pickle                                                           
import scipy.io as sio                                                  
from load_objects import *                                              
from mpl_toolkits.basemap import Basemap
import numpy as np

                                                                        
def plot_berg_location(xil,yil,show=False):
    """
    This function plots the xy position of each iceberg over time.

    Args:
        xil (list): X location of the iceberg.
        yil (list): Y location of the iceberg.

    Generates:
        A plot showing the xy position of icebergs over time.
    """
    f = plt.figure()
    plt.plot(xil.transpose(),yil.transpose())                                                       
    plt.xlabel('X Location')                                                   
    plt.ylabel('Y Location')                                                   
    plt.title('Iceberg Drift')                                          
    if show:
        plt.show()
    else:
        return f        

def plot_track_on_map(berg_lon,berg_lat,show=False):

    berg_lon = berg_lon[~np.isnan(berg_lon)]-360
    berg_lat = berg_lat[~np.isnan(berg_lat)]

    # Lambert Conformal Conic map.

    llc_lon = np.amin(berg_lon) - 10
    llc_lat = np.amin(berg_lat) - 10

    urc_lon = np.amax(berg_lon) + 10
    urc_lat = np.amax(berg_lat) + 10

    fig = plt.figure()
    
    m = Basemap(llcrnrlon=llc_lon,llcrnrlat=llc_lat,urcrnrlon=urc_lon,urcrnrlat=urc_lat,
                projection='lcc',lat_1=20., lon_0=-60.,
                resolution ='l',area_thresh=1000.)

    x, y = m(berg_lon,berg_lat)
    m.plot(x,y,linewidth=1.5,color='r')

    # draw coastlines, meridians and parallels.
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    m.drawparallels(np.arange(10,70,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(-100,0,20),labels=[0,0,0,1])
    plt.title('Iceberg Drift')
    if show:
        plt.show()
    else:
        return fig
