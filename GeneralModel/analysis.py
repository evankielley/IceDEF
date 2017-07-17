"""
This module does analyis on the output of an Iceberg simulation.
"""

import matplotlib.pyplot as plt                                         
import pickle                                                           
import scipy.io as sio                                                  
from load_objects import *                                              
                                                                        
def compare_outputs(mXIL,mYIL,bb,pyOutloc,trajnum,nt,relative=False):                  
    """
    This function compares the output of two separate models.

    Args: 
        mXIL (list): Matlab iceberg x location matrix
        mYIL (list): Matlab iceberg y location matrix
        bb (int): Iceberg size class
        pyOutloc (str): Path to Python model output
        trajnum (int): Number of trajectories per iceberg size class
        nt (int): Number of timesteps per run


    Returns:
        dXIL (list): Difference between X locations between model outputs.
        dYIL (list): Difference between Y locations between model outputs.
    """
    pyXIL, pyYIL = load_objects(pyOutloc,trajnum,nt)                    
    dXIL = np.empty([trajnum,nt])*np.nan    
    dYIL = np.empty([trajnum,nt])*np.nan    

    if relative:
        for j in range(0,nt):
            #xlen = np.where(np.isnan(pyXIL[i,:]))[0][-1]
            #ylen = np.where(np.isnan(pyYIL[i,:]))[0][-1]
            dXIL[:,j] = np.divide(mXIL[bb-1,:,j]-pyXIL[:,j],j)                                  
            dYIL[:,j] = np.divide(mYIL[bb-1,:,j]-pyYIL[:,j],j)                                  
        return dXIL, dYIL

    else:
        dXIL = mXIL[bb-1,:,:] - pyXIL[:,:]                                  
        dYIL = mYIL[bb-1,:,:] - pyYIL[:,:]                                  
        return dXIL, dYIL 

def plot_model_diff(xDiff,yDiff,show=False):
    """
    This function plot the X and Y location difference between models.

    Args:
        xDiff (list): A list of all the X location differences between models.
        yDiff (list): A list od all the Y location differences between models.

    Generates:
        A plot with two subplots; one for x and y differences; respectively.
    """
    f = plt.figure()
    plt.subplot(121)
    plt.plot(xDiff.transpose())
    plt.ylabel('Difference in Location')
    plt.xlabel('Timestep')
    plt.title('X')
    plt.subplot(122)
    plt.plot(yDiff.transpose())
    #plt.ylabel('Difference in Y Location')
    plt.xlabel('Timestep')
    plt.title('Y')
    if show:
        plt.show()
    else:
        return f        

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

from mpl_toolkits.basemap import Basemap

def plot_track_on_map(xil,yil,show=False):
    xil = xil - 360
    xil = xil[0]
    xil = xil[:360]
    print(xil)
    yil = yil[0]
    yil = yil[:360]
    print(yil)
    # Lambert Conformal Conic map.
    m = Basemap(llcrnrlon=-80.,llcrnrlat=30.,urcrnrlon=-30.,urcrnrlat=70.,
                projection='lcc',lat_1=20., lon_0=-60.,
                resolution ='l',area_thresh=1000.)
    #m = Basemap(llcrnrlon=-80.,llcrnrlat=30.,urcrnrlon=-30.,urcrnrlat=80.,
    #            projection='lcc',lat_1=20.,lat_2=40.,lon_0=-60.,
    #            resolution ='l',area_thresh=1000.)

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

def plot_track_on_map2(xil,yil,show=False):
    xil = xil - 360
    xil = xil[0]
    xil = xil[:30]
    print(xil)
    yil = yil[0]
    yil = yil[:30]
    print(yil)
    # draw map with markers for float locations
    m = Basemap(projection='hammer',lon_0=0)
    x, y = m(xil,yil)
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    m.scatter(x,y,3,marker='o',color='k')
    plt.title('Iceberg Drift')
    plt.show()



