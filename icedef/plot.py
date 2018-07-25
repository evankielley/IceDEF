"""This module makes plots and animations for objects in the icedef package.
"""

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.animation import FuncAnimation
import numpy as np
import netCDF4 as nc
from scipy.interpolate import RegularGridInterpolator as RGI
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


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

def get_mercator_basemap(min_lon, max_lon, min_lat, max_lat):

    num_ticks = 5
    lat_pad = abs(min_lat - max_lat)/num_ticks
    lon_pad = abs(min_lon - max_lon)/num_ticks

    # instantiate Basemap object with mercator projection
    m = Basemap(projection='merc', lat_0 = 57, lon_0 = -135,
                resolution = 'l', area_thresh = 0.1,
                llcrnrlon = min_lon - lon_pad,
                llcrnrlat = min_lat - lat_pad,
                urcrnrlon = max_lon + lon_pad,
                urcrnrlat = max_lat + lat_pad)


    # parallels are lines of latitude, meridians are lines of longitude
    parallels = np.round(np.arange(min_lat, max_lat + lat_pad, lat_pad), 2)
    meridians = np.round(np.arange(min_lon, max_lon + lon_pad, lon_pad), 2)

    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[True,False,False,False])
    mers = m.drawmeridians(meridians,labels=[False,False,False,True])
    for mer in mers:
        try:
            mers[mer][1][0].set_rotation(90)
        except:
            pass

    # make the land pretty
    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color = 'coral')
    m.drawmapboundary()

    return m

def plot_drift_track_test_case(iip_berg, mod_berg, time_labels=True):
    """This function plots the drift track of a simulated iceberg against its observed coordinates"

    Args:
        iip_berg (Iceberg): all attributes correspond to data from a particular IIP iceberg.
        mod_berg (Iceberg): initial attributes and final time correspond to data from a particular
                            IIP iceberg but the rest comes from drift simulation.
    """

    fig, ax = plt.subplots(dpi=150)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Comparison of Drift Model Track to IIP Observations')
    ax.xaxis.labelpad = 60
    ax.yaxis.labelpad = 60

    iip_lons = iip_berg.history['X']
    iip_lats = iip_berg.history['Y']
    mod_lons = mod_berg.history['X']
    mod_lats = mod_berg.history['Y']

    # Get min and max iceberg lons and lats for defining the plot axes ranges
    berg_xmin = min(min(mod_lons), min(iip_lons))
    berg_xmax = max(max(mod_lons), max(iip_lons))
    berg_ymin = min(min(mod_lats), min(iip_lats))
    berg_ymax = max(max(mod_lats), max(iip_lats))

    # spatial buffer for expanding the plot axes ranges for better readability
    buff = 0.1

    m = get_mercator_basemap(berg_xmin, berg_xmax, berg_ymin, berg_ymax)

    # map lons and lats to x and y in meters
    iip_x, iip_y = m(iip_lons, iip_lats)
    mod_x, mod_y = m(mod_lons, mod_lats)

    # scatter both sets of data points so timesteps can be seen
    ax.scatter(iip_x, iip_y, marker='o', s=10, c='black')
    ax.scatter(mod_x, mod_y, marker='o', s=1, c='blue')

    # if True, annotate matching points in time between datasets
    if time_labels:

        iip_times = iip_berg.history['T']
        mod_times = mod_berg.history['T']

        matching_indices = find_matching_value_indices(iip_times, mod_times)

        iip_hours = [(t-iip_times[0]).days*24 + (t-iip_times[0]).seconds/3600 for t in iip_times]
        hour_labels = [str(round(x, 1)) for x in iip_hours]

        for i, hour_label in enumerate(hour_labels):
            ax.text(iip_x[i], iip_y[i], hour_label)

        for i, hour_label in enumerate(hour_labels):

            try:
                j = matching_indices[i]
            except IndexError:
                break

            if not round(iip_lons[i], 2) == round(mod_lons[j], 2):
                ax.scatter(mod_x[j], mod_y[j], marker='o', color='red')
                ax.text(mod_x[j], mod_y[j], hour_label)

    return fig



def plot_spaghetti_test_case(iip_berg, mod_berg_list, time_labels=False):

    fig, ax = plt.subplots(dpi=150)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Comparison of Drift Model Tracks to IIP Observations')
    ax.xaxis.labelpad = 60
    ax.yaxis.labelpad = 60

    iip_lons = iip_berg.history['X']
    iip_lats = iip_berg.history['Y']

    berg_xmin, berg_xmax = min(iip_lons), max(iip_lons)
    berg_ymin, berg_ymax = min(iip_lats), max(iip_lats)

    for mod_berg in mod_berg_list:

        if not mod_berg.out_of_bounds:

            mod_lons = mod_berg.history['X']
            mod_lats = mod_berg.history['Y']

            # Get min and max iceberg lons and lats for defining the plot axes ranges
            tmp_xmin = min(min(mod_lons), min(iip_lons))
            tmp_xmax = max(max(mod_lons), max(iip_lons))
            tmp_ymin = min(min(mod_lats), min(iip_lats))
            tmp_ymax = max(max(mod_lats), max(iip_lats))

            berg_xmin = min(berg_xmin, tmp_xmin)
            berg_xmax = max(berg_xmax, tmp_xmax)
            berg_ymin = min(berg_ymin, tmp_ymin)
            berg_ymax = max(berg_ymax, tmp_ymax)


    m = get_mercator_basemap(berg_xmin, berg_xmax, berg_ymin, berg_ymax)


    # map lons and lats to x and y in meters
    iip_x, iip_y = m(iip_lons, iip_lats)

    # get lists of iip times for comparing to mod times to later annotate plot with labels
    iip_times = iip_berg.history['T']
    iip_hours = [(t-iip_times[0]).days*24 + (t-iip_times[0]).seconds/3600 for t in iip_times]
    hour_labels = [str(round(x, 1)) for x in iip_hours]


    for i, hour_label in enumerate(hour_labels):
        ax.scatter(iip_x[i], iip_y[i], marker='o', color='black')
        ax.text(iip_x[i], iip_y[i], hour_label)


    for mod_berg in mod_berg_list:

        mod_lons = mod_berg.history['X']
        mod_lats = mod_berg.history['Y']
        mod_x, mod_y = m(mod_lons, mod_lats)
        ax.scatter(mod_x, mod_y, marker='o', s=1)#, c='red')

        mod_times = mod_berg.history['T']

        matching_indices = find_matching_value_indices(iip_times, mod_times)

        for i, hour_label in enumerate(hour_labels):

            try:
                j = matching_indices[i]
            except IndexError:
                break

            if not round(iip_lons[i], 2) == round(mod_lons[j], 2):
                ax.scatter(mod_x[j], mod_y[j], marker='o', color='black',s=10)
                if time_labels:
                    ax.text(mod_x[j], mod_y[j], hour_label)


    plt.savefig('test_plot.png')

    return fig


def berg_metocean_animation(xyt_berg, xyt_grid, grid_scales, uv_data,
                             fname='field_anim', v_auto = True, vmin=0, vmax=1,
                             scale=1, headwidth=5, width=5e-3, speed=100):

    x_berg, y_berg, t_berg = xyt_berg
    x_grid, y_grid, t_grid = xyt_grid
    dx_scale, dy_scale, dt_scale = grid_scales
    u_data, v_data = uv_data

    dx = np.mean(np.diff(x_grid))
    dy = np.mean(np.diff(y_grid))
    dt = np.mean(np.diff(t_grid))

    x0 = min(x_berg) - dx
    xf = max(x_berg) + dx
    y0 = min(y_berg) - dy
    yf = max(y_berg) + dy
    t0 = min(t_berg) - dt
    tf = max(t_berg) + dt

    xid0 = np.where(abs(x_grid - x0) < dx)[0][0]
    xidf = np.where(abs(x_grid - xf) < dx)[0][-1]
    yid0 = np.where(abs(y_grid - y0) < dy)[0][0]
    yidf = np.where(abs(y_grid - yf) < dy)[0][-1]
    tid0 = np.where(abs(t_grid - t0) < dt)[0][0]
    tidf = np.where(abs(t_grid - tf) < dt)[0][-1]

    # overwrite existing variables with a subset of them
    x_grid = x_grid[xid0:xidf+1]
    y_grid = y_grid[yid0:yidf+1]
    t_grid = t_grid[tid0:tidf]
    u_data = u_data[tid0:tidf, yid0:yidf+1, xid0:xidf+1]
    v_data = v_data[tid0:tidf, yid0:yidf+1, xid0:xidf+1]

    # make interpolators for u and v data
    u_RGI = RGI((t_grid, y_grid, x_grid), u_data)
    v_RGI = RGI((t_grid, y_grid, x_grid), v_data)

    # overwrite again with scales
    x_grid = np.arange(x_grid[0], x_grid[-1], dx/dx_scale)
    y_grid = np.arange(y_grid[0], y_grid[-1], dy/dy_scale)
    t_grid = np.arange(t_grid[0], t_grid[-1], dt/dt_scale)

    # initialize matrices
    u_mat = np.empty([len(t_grid), len(y_grid), len(x_grid)])
    v_mat = np.empty([len(t_grid), len(y_grid), len(x_grid)])

    # fill matrices with interpolated values
    for i, ival in enumerate(t_grid):
        for j, jval in enumerate(y_grid):
            for k, kval in enumerate(x_grid):
                u_mat[i][j][k] = u_RGI([ival,jval,kval])
                v_mat[i][j][k] = v_RGI([ival,jval,kval])

    # make matrix of magnitudes of u and v
    w_mat = np.sqrt(u_mat**2 + v_mat**2)

    # set colorbar bounds
    if v_auto:
        vmin = w_mat.min()
        vmax = w_mat.max()

    # setup figure, background, axes, and gridlines
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([x0, xf, y0, yf], ccrs.PlateCarree())
    ax.text(-0.2, 0.5, 'latitude', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=ax.transAxes)
    ax.text(0.5, -0.2, 'longitude', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',
            transform=ax.transAxes)
    ax.stock_img()
    ax.coastlines('50m')
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=2,color='gray',alpha=0.5,linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # plot
    im = plt.imshow(w_mat[0,:,:], extent=[x0, xf, y0, yf], origin = 'lower', vmin=vmin, vmax=vmax)
    line, = plt.plot(x_berg[0], y_berg[0], color='black')
    plt.colorbar()  # colorbar must be called before quiv
    quiv = plt.quiver(x_grid, y_grid, u_mat[0,:,:], v_mat[0,:,:], scale=scale, headwidth=headwidth, width=width)
    title = plt.title('time: 0 hours')

    def animate(i):
        im.set_data(w_mat[i,:,:])
        quiv.set_UVC(u_mat[i,:,:], v_mat[i,:,:])
        title.set_text('time: {:.0f} hours'.format(t_grid[i]-t_grid[0]))
        min_indx = np.argmin([abs(t) for t in t_berg - t_grid[i]])
        line.set_data(x_berg[0:min_indx+1], y_berg[0:min_indx+1])
        return im, line

    anim = FuncAnimation(fig, animate, frames=w_mat[:,0,0].size-1, interval=speed)
    #HTML(anim.to_html5_video())
    anim.save(f'plots/{fname}.gif',writer='imagemagick')
