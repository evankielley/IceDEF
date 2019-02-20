"""This module makes plots and animations for objects in the icedef package.
"""

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.animation import FuncAnimation
import numpy as np
import xarray as xr
from pandas import to_datetime

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)   # fontsize of the figure title


def get_map(**kwargs):
    """This function returns a Basemap map.

    Args:
        **kwargs: coming soon.

    Returns:
        map_ (mpl_toolkits.basemap.Basemap): a Basemap map.

    """

    projection = kwargs.pop('projection', 'stere')
    lon_0 = kwargs.pop('lon_0', -50)
    lat_0 = kwargs.pop('lat_0', 50)
    lat_ts = kwargs.pop('lat_ts', 45)
    resolution = kwargs.pop('resolution', 'l')
    area_thresh = kwargs.pop('area_thresh', 0.1)
    llcrnrlon = kwargs.pop('llcrnrlon', -60)
    llcrnrlat = kwargs.pop('llcrnrlat', 40)
    urcrnrlon = kwargs.pop('urcrnrlon', -40)
    urcrnrlat = kwargs.pop('urcrnrlat', 60)

    map_ = Basemap(projection=projection,
                   lon_0=lon_0, lat_0=lat_0, lat_ts=lat_ts,
                   resolution=resolution, area_thresh=area_thresh,
                   llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)

    return map_


def draw_map(map_, **kwargs):
    """This function draws background information like labels, coastlines, and grid lines on a Basemap map.

    Args:
        map_ (mpl_toolkits.basemap.Basemap): a Basemap map.
        **kwargs: coming soon.

    """

    drawcoastlines = kwargs.pop('drawcoastlines', True)
    drawstates = kwargs.pop('drawstates', False)
    drawcountries = kwargs.pop('drawcountries', False)
    parallels = kwargs.pop('parallels', np.arange(0, 90, 5))
    meridians = kwargs.pop('meridians', np.arange(0., 360., 10.))
    xtick_rotation_angle = kwargs.pop('xtick_rotation_angle', 90)

    if drawcoastlines:
        map_.drawcoastlines()
    if drawstates:
        map_.drawstates()
    if drawcountries:
        map_.drawcountries()

    # draw parallels
    map_.drawparallels(parallels, labels=[1, 0, 0, 0])

    # draw meridians
    meridians = map_.drawmeridians(meridians, labels=[0, 0, 0, 1])

    if xtick_rotation_angle:

        # rotate x tick labels
        for meridian in meridians:

            try:
                meridians[meridian][1][0].set_rotation(xtick_rotation_angle)

            except IndexError:
                pass


def plot_track(*latlons, **kwargs):

    min_lat = None
    min_lon = None
    max_lat = None
    max_lon = None

    for latlon in latlons:

        lats, lons = latlon
        temp_min_lat = min(lats)
        temp_min_lon = min(lons)
        temp_max_lat = max(lats)
        temp_max_lon = max(lons)
        if min_lat is None or temp_min_lat < min_lat:
            min_lat = temp_min_lat
        if min_lon is None or temp_min_lon < min_lon:
            min_lon = temp_min_lon
        if max_lat is None or temp_max_lat > max_lat:
            max_lat = temp_max_lat
        if max_lon is None or temp_max_lon > max_lon:
            max_lon = temp_max_lon

    scatter_kwargs = kwargs.pop('scatter_kwargs', {})
    quiver_kwargs = kwargs.pop('quiver_kwargs', {})
    map_kwargs = kwargs.pop('map_kwargs', {})
    legend_kwargs = kwargs.pop('legend_kwargs', {})

    new_map_kwargs = get_map_kwargs(min_lat, min_lon, max_lat, max_lon, **map_kwargs)
    map_kwargs.update(new_map_kwargs)

    fig, ax = plt.subplots()
    map_ = get_map(**map_kwargs)
    draw_map(map_, **map_kwargs)

    labels = kwargs.pop('labels', None)
    quiver_labels = quiver_kwargs.get('label', None)

    if labels is not None:
        show_legend = True
        quiver_kwargs['label'] = [''] * 100

    elif quiver_labels is not None:
        show_legend = True
        labels = [''] * 100

    else:
        show_legend = False
        labels = [''] * 100
        quiver_kwargs['label'] = [''] * 100

    i = 0

    for latlon in latlons:

        lats, lons = latlon
        eastings, northings = map_(lons, lats)
        ax.scatter(eastings, northings, label=labels[i], **scatter_kwargs)
        i += 1

    vectors = kwargs.pop('vectors', None)

    if vectors is not None:

        ax = plot_quivers(eastings, northings, vectors, ax, **quiver_kwargs)

    title = kwargs.pop('title', '')
    ax.set_title(title)
    if show_legend:
        ax.legend(**legend_kwargs)

    return fig, ax


def plot_quivers(x, y, vectors, ax=None, **kwargs):

    gap = kwargs.pop('gap', 10)
    color = kwargs.pop('color', ['black'] * 100)
    label = kwargs.pop('label', [''] * 100)

    if ax is None:

        ax = plt.subplot(111)

    i = 0

    for vector_u, vector_v in vectors:

        if isinstance(vector_u, xr.core.dataarray.DataArray):
            vector_u = vector_u.values

        if isinstance(vector_v, xr.core.dataarray.DataArray):
            vector_v = vector_v.values

        ax.quiver(x[::gap], y[::gap], vector_u[::gap], vector_v[::gap], color=color[i], label=label[i], **kwargs)

        i += 1

    return ax


def get_map_kwargs(min_lat, min_lon, max_lat, max_lon, **kwargs):

    pads = kwargs.pop('pads', [0.01] * 4)

    min_lat_padded = min_lat - pads[1]
    max_lat_padded = max_lat + pads[3]
    min_lon_padded = min_lon - pads[2]
    max_lon_padded = max_lon + pads[0]

    map_kwargs = {

        'projection': kwargs.pop('projection', 'stere'),
        'llcrnrlon': kwargs.pop('llcrnrlon', min_lon_padded),
        'urcrnrlon': kwargs.pop('urcrnrlon', max_lon_padded),
        'llcrnrlat': kwargs.pop('llcrnrlat', min_lat_padded),
        'urcrnrlat': kwargs.pop('urcrnrlat', max_lat_padded),
        'parallels': kwargs.pop('parallels', np.linspace(min_lat_padded, max_lat_padded + pads[1], 10)),
        'meridians': kwargs.pop('meridians', np.linspace(min_lon_padded, max_lon_padded + pads[0], 10)),
        'lon_0': kwargs.pop('lon_0', (min_lon + max_lon) / 2),
        'lat_0': kwargs.pop('lat_0', (min_lat + max_lat) / 2),
        'lat_ts': kwargs.pop('lat_ts', (min_lat + max_lat) / 2),
    }

    return map_kwargs


def plot_image(lats, lons, data, **kwargs):

    min_lat = min(lats)
    min_lon = min(lons)
    max_lat = max(lats)
    max_lon = max(lons)

    map_kwargs = get_map_kwargs(min_lat, min_lon, max_lat, max_lon, **kwargs)

    kwargs.update(map_kwargs)

    map_ = get_map(**kwargs)

    pcolormesh = True
    imshow = False

    if pcolormesh:

        meshed_lons, meshed_lats = np.meshgrid(lons, lats)
        x, y = map_(meshed_lons, meshed_lats)
        map_.pcolormesh(x, y, data)

    elif imshow:

        x0, y0 = map_(min_lon, min_lat)
        x1, y1 = map_(max_lon, max_lat)
        map_.imshow(data, extent=[x0, x1, y0, y1], origin='lower')


def animate_field(lats, lons, times, u, v, gap=10, time_format='%Y-%m-%d %H'):

    def dt64_to_str(t): return to_datetime(t).strftime(time_format)

    fig = plt.figure(figsize=(10, 10))
    map_ = get_map()

    magnitude = np.sqrt(u ** 2 + v ** 2)

    image = map_.imshow(magnitude[0, :, :], origin='lower')
    plt.colorbar()
    title = plt.title(dt64_to_str(times[0]))

    meshed_lons, meshed_lats = np.meshgrid(lons[::gap], lats[::gap])
    x, y = map_(meshed_lons, meshed_lats)
    quiver = map_.quiver(x, y, u[0, ::gap, ::gap], v[0, ::gap, ::gap])

    def animate(i):
        image.set_data(magnitude[i, :, :])
        title.set_text(dt64_to_str(times[i]))
        quiver.set_UVC(u[i, ::gap, ::gap], v[i, ::gap, ::gap])

    return FuncAnimation(fig, animate, frames=len(times))
