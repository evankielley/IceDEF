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


def plot_iceberg_track(lats, lons, **kwargs):
    plot_width = kwargs.pop('plot_width', 10)
    plot_height = kwargs.pop('plot_height', 8)

    fig = plt.figure(figsize=(plot_width, plot_height))

    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    vectors = kwargs.pop('vectors', None)
    gap = kwargs.pop('gap', 10)
    track_marker_size = kwargs.pop('track_marker_size', 5)
    ref_marker_size = kwargs.pop('ref_marker_size', 5)
    arrow_scale = kwargs.pop('arrow_scale', None)
    arrow_shaftwidth = kwargs.pop('arrow_shaftwidth', 0.0005 * plot_width)
    arrow_headlength = kwargs.pop('arrow_headlength', 5)
    arrow_headwidth = kwargs.pop('arrow_headwidth', 3)
    arrow_colors = kwargs.pop('arrow_colors', ['black'] * 10)
    arrow_labels = kwargs.pop('arrow_labels', [''] * 10)
    ref_track = kwargs.pop('ref_track', None)
    quiver_ref_track = kwargs.pop('quiver_ref_track', False)

    plt.scatter(lons, lats, s=track_marker_size, label='')

    if vectors is not None:

        i = 0

        for vector_u, vector_v in vectors:
            plt.quiver(lons[::gap], lats[::gap], vector_u[::gap], vector_v[::gap],
                       scale=arrow_scale, width=arrow_shaftwidth, headlength=arrow_headlength,
                       headwidth=arrow_headwidth, color=arrow_colors[i], label=arrow_labels[i])
            i += 1

        if ref_track is not None:

            ref_lats, ref_lons = ref_track
            plt.scatter(ref_lons, ref_lats, s=ref_marker_size, label='')
            plt.plot(ref_lons, ref_lats)

            if quiver_ref_track:

                ref_lat_list = [];
                ref_lon_list = []

                for i in np.arange(0, len(vector_u), gap):
                    ref_lat_list.append(ref_lats.interp(time=xds.time[i]))
                    ref_lon_list.append(ref_lons.interp(time=xds.time[i]))

                i = 0

                for vector_u, vector_v in vectors:
                    plt.quiver(ref_lon_list, ref_lat_list, vector_u[::gap], vector_v[::gap],
                               scale=arrow_scale, width=arrow_shaftwidth, headlength=arrow_headlength,
                               headwidth=arrow_headwidth, color=arrow_colors[i], label='')
                    i += 1

        if arrow_labels[0] is not '':
            plt.legend()

    plt.show()


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
    m.fillcontinents(color='coral')
    m.drawmapboundary()

    return m


def iceberg_metocean_animation(xyt_berg, xyt_grid, grid_scales, uv_data,
                             fname='field_anim', v_auto=True, vmin=0, vmax=1,
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
