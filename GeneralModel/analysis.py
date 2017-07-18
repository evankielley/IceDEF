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

def plot_track_on_map(xil,yil,show=False):

    #xil = xil - 360
    xil = xil - 360
    xil = xil[~np.isnan(xil)]
    #yil = yil[0]
    yil = yil[~np.isnan(yil)]

    # Lambert Conformal Conic map.
    m = Basemap(llcrnrlon=-80.,llcrnrlat=30.,urcrnrlon=-30.,urcrnrlat=70.,
                projection='lcc',lat_1=20., lon_0=-60.,
                resolution ='l',area_thresh=1000.)

    x, y = m(xil,yil)
    m.plot(x,y,linewidth=1.5,color='r')

    # draw coastlines, meridians and parallels.
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    m.drawparallels(np.arange(10,70,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(-100,0,20),labels=[0,0,0,1])
    plt.title('Iceberg Drift')
    plt.show()
