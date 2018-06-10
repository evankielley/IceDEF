"""This module makes plots and animations for objects in the icedef package.
"""

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.animation import FuncAnimation
import numpy as np
import netCDF4 as nc


def find_matching_value_indices(small_list, big_list):
    """This function finds the indices from one list of the values that match the values from the other list.
    
    Args:
        small_list (list of float or int): 1D list of numbers, smaller in length than big_list
        big_list (list of float or int): 1D list of evenly spaced numbers, larger in length than small_list
        
    Returns:
        matching_indices (list of int): list of all the indices of the bigger list which match the smaller list values
        
    """
    
    big_delta = big_list[1] - big_list[0]
    
    matching_indices = []
    
    for small_item in small_list:
        for big_index, big_item in enumerate(big_list):
            delta = abs(small_item - big_item)
            if delta < big_delta:
                matching_index = big_index
                matching_indices.append(matching_index)
                break
    
    return matching_indices



def plot_drift_track_test_case(iip_berg, mod_berg):
    """This function plots the drift track of a simulated iceberg against its observed coordinates"
    
    Args:
        iip_berg (Iceberg): all attributes correspond to data from a particular IIP iceberg.
        mod_berg (Iceberg): initial attributes and final time correspond to data from a particular 
                            IIP iceberg but the rest comes from drift simulation.  
    """
    
    buff = 0.1

    iip_lons = iip_berg.history['X']
    iip_lats = iip_berg.history['Y']
    mod_lons = mod_berg.history['X']
    mod_lats = mod_berg.history['Y']
    
    berg_xmin = min(min(mod_lons), min(iip_lons)) 
    berg_xmax = max(max(mod_lons), max(iip_lons))
    berg_ymin = min(min(mod_lats), min(iip_lats))
    berg_ymax = max(max(mod_lats), max(iip_lats))

    m = Basemap(projection='merc', lat_0 = 57, lon_0 = -135,
                resolution = 'l', area_thresh = 0.1,
                llcrnrlon = berg_xmin - buff, 
                llcrnrlat = berg_ymin - buff,
                urcrnrlon = berg_xmax + buff, 
                urcrnrlat = berg_ymax + buff)
    
    # parallels are lines of latitude, meridians are lines of longitude (labels = [left,right,top,bottom])
    parallels = np.arange(berg_ymin, berg_ymax + buff, buff)
    m.drawparallels(parallels,labels=[True,False,False,False])
    meridians = np.arange(berg_xmin, berg_xmax + buff, buff)
    m.drawmeridians(meridians,labels=[False,False,False,True])

    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color = 'coral')
    m.drawmapboundary()

    iip_x, iip_y = m(iip_lons, iip_lats)
    mod_x, mod_y = m(mod_lons, mod_lats)
        
    m.plot(iip_x, iip_y, color='black') 
    m.plot(mod_x, mod_y, 'bo', markersize=1)
    
    iip_times = iip_berg.history['T']
    mod_times = mod_berg.history['T']

    matching_indices = find_matching_value_indices(iip_times, mod_times)
                
    for i in matching_indices:
        tdelta = mod_times[i] - mod_times[0]
        hour_label = tdelta.days*24 + tdelta.seconds/3600
        plt.scatter(mod_x[i], mod_y[i], color='r')
        plt.text(mod_x[i], mod_y[i], hour_label)
    
    iip_hours = [(t-iip_times[0]).days*24 + (t-iip_times[0]).seconds/3600 for t in iip_times]
    plt.scatter(iip_x, iip_y, color='r')
    plt.text(iip_x, mod_y, iip_hours)

    plt.show()


        
def animate_currents(ocean_data, iip_berg, mod_berg):
    """This function saves a gif of the ocean current velocities alongside an iceberg drift test case.
    
    Args:
        ocean_data (ECMWF_Ocean): ocean data for time and space of iceberg drift simulation.
        iip_berg (Iceberg): all attributes correspond to data from a particular IIP iceberg.
        mod_berg (Iceberg): initial attributes and final time correspond to data from a particular 
                            IIP iceberg but the rest comes from drift simulation.
    """
    
    
    mod_berg_t1950 = nc.date2num(mod_berg.datetimes, 
                                 'hours since 1950-01-01 00:00:00.0 00:00', 'standard')
    
    lon0 = np.where(ocean_data.lons < min(mod_berg.lons))[0][-1]
    lonn = np.where(ocean_data.lons > max(mod_berg.lons))[0][0]
    lat0 = np.where(ocean_data.lats < min(mod_berg.lats))[0][-1] 
    latn = np.where(ocean_data.lats > max(mod_berg.lats))[0][0]
    
    UW = np.empty([len(mod_berg_t1950), 
                len(ocean_data.lats[lat0-1:latn+1]), len(ocean_data.lons[lon0-1:lonn+1])])
    VW = np.empty([len(mod_berg_t1950), 
                len(ocean_data.lats[lat0-1:latn+1]), len(ocean_data.lons[lon0-1:lonn+1])])
    
    for i, ival in enumerate(mod_berg_t1950):
        for j, jval in enumerate(ocean_data.lats[lat0-1:latn+1]):
            for k, kval in enumerate(ocean_data.lons[lon0-1:lonn+1]):
                UW[i][j][k] = ocean_data.iUW([ival,jval,kval])
                VW[i][j][k] = ocean_data.iVW([ival,jval,kval])
    
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([ocean_data.lons[lon0], ocean_data.lons[lonn], ocean_data.lats[lat0], ocean_data.lats[latn]], ccrs.PlateCarree())
    #ax.stock_img()
    ax.coastlines('50m')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ax.scatter(iip_berg.lons[:], iip_berg.lats[:], color='black')
    #ax.plot(mod_berg.lons[:], mod_berg.lats[:], color='black')
    water_mag = np.sqrt(UW**2 + VW**2)
    im = plt.imshow(water_mag[0,:,:], 
                    extent=[ocean_data.lons[lon0-1], ocean_data.lons[lonn+1], 
                            ocean_data.lats[lat0-1], ocean_data.lats[latn+1]],
                    origin = 'lower', vmin=0, vmax=0.5)
    
    line, = plt.plot(mod_berg.lons[0], mod_berg.lats[0], color='black')
    
    plt.colorbar()
    quiv = plt.quiver(ocean_data.lons[lon0-1:lonn+1], ocean_data.lats[lat0-1:latn+1], UW[0,:,:], VW[0,:,:], 
                      scale=1, headwidth=5, width=0.005)
    title = plt.title('time: 0 hours')


    def animate(i):

        im.set_data(water_mag[i,:,:])
        quiv.set_UVC(UW[i,:,:], VW[i,:,:])
        title.set_text('time: {:.0f} hours'.format(mod_berg_t1950[i]-mod_berg_t1950[0]))
        line.set_data(mod_berg.lons[0:i+1], mod_berg.lats[0:i+1])
        
        
        return im, line
    
    anim = FuncAnimation(fig, animate, frames=water_mag[:,0,0].size-1, interval=100)
    #HTML(anim.to_html5_video())
    anim.save('plots/water_mag.gif',writer='imagemagick')
        
    
    
def animate_winds(atm_data, iip_berg, mod_berg):
    """This function saves a gif of the wind velocities alongside an iceberg drift test case.
    
    Args:
        atm_data (ECMWF_Atm): atmospheric data for time and space of iceberg drift simulation.
        iip_berg (Iceberg): all attributes correspond to data from a particular IIP iceberg.
        mod_berg (Iceberg): initial attributes and final time correspond to data from a particular 
                            IIP iceberg but the rest comes from drift simulation.
    """
    
    mod_berg_t1900 = nc.date2num(mod_berg.datetimes, 
                                 'hours since 1900-01-01 00:00:00.0 00:00', 'standard')
    
    lon0 = np.where(atm_data.lons < min(mod_berg.lons))[0][-1]
    lonn = np.where(atm_data.lons > max(mod_berg.lons))[0][0]
    lat0 = np.where(atm_data.lats < min(mod_berg.lats))[0][-1] 
    latn = np.where(atm_data.lats > max(mod_berg.lats))[0][0]
    
    UA = np.empty([len(mod_berg_t1900), 
                len(atm_data.lats[lat0-1:latn+1]), len(atm_data.lons[lon0-1:lonn+1])])
    VA = np.empty([len(mod_berg_t1900), 
                len(atm_data.lats[lat0-1:latn+1]), len(atm_data.lons[lon0-1:lonn+1])])
    
    for i, ival in enumerate(mod_berg_t1900):
        for j, jval in enumerate(atm_data.lats[lat0-1:latn+1]):
            for k, kval in enumerate(atm_data.lons[lon0-1:lonn+1]):
                UA[i][j][k] = atm_data.iUA([ival,jval,kval])
                VA[i][j][k] = atm_data.iVA([ival,jval,kval])
    
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([atm_data.lons[lon0], atm_data.lons[lonn], atm_data.lats[lat0], atm_data.lats[latn]], ccrs.PlateCarree())
    #ax.stock_img()
    ax.coastlines('50m')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ax.scatter(iip_berg.lons[:], iip_berg.lats[:], color='black')
    wind_mag = np.sqrt(UA**2 + VA**2)
    im = plt.imshow(wind_mag[0,:,:], 
                    extent=[atm_data.lons[lon0-1], atm_data.lons[lonn+1], 
                            atm_data.lats[lat0-1], atm_data.lats[latn+1]],
                    origin = 'lower', vmin=5, vmax=15)
    
    line, = plt.plot(mod_berg.lons[0], mod_berg.lats[0], color='black')
    
    plt.colorbar()
    quiv = plt.quiver(atm_data.lons[lon0-1:lonn+1], atm_data.lats[lat0-1:latn+1], UA[0,:,:], VA[0,:,:], 
                      scale=20, headwidth=5, width=0.005)
    title = plt.title('time: 0 hours')


    def animate(i):

        im.set_data(wind_mag[i,:,:])
        quiv.set_UVC(UA[i,:,:], VA[i,:,:])
        title.set_text('time: {:.0f} hours'.format(mod_berg_t1900[i]-mod_berg_t1900[0]))
        line.set_data(mod_berg.lons[0:i+1], mod_berg.lats[0:i+1])
        
        
        return im, line
    
    anim = FuncAnimation(fig, animate, frames=wind_mag[:,0,0].size-1, interval=100)
    #HTML(anim.to_html5_video())
    anim.save('plots/wind_mag.gif',writer='imagemagick')