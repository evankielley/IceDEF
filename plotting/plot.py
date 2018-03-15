import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.basemap import Basemap
import bisect
from matplotlib.animation import FuncAnimation


def plot0(iip_berg, mod_berg):
    tol = 0.1  # tstep
    for i,t in enumerate(iip_berg.t2000):
        diff = abs(np.asarray(mod_berg.t2000) - t)
        for j,diff_t in enumerate(diff):
            if diff_t < tol:
                break
        plt.plot(iip_berg.lats[i], iip_berg.lons[i], marker='o', color='red')
        plt.text(iip_berg.lats[i], iip_berg.lons[i], '{0:.1f}'.format(t-iip_berg.t2000[0]))
        plt.plot(mod_berg.lats[j], mod_berg.lons[j], marker='o', color='yellow')
        plt.text(mod_berg.lats[j], mod_berg.lons[j], '{0:.1f}'.format(t-iip_berg.t2000[0]))
    plt.plot(iip_berg.lats, iip_berg.lons, label='observed', color='red')
    plt.plot(mod_berg.lats, mod_berg.lons, label='computed', color='orange')
    plt.legend()
    plt.xlabel('Latitude'); plt.ylabel('Longitude')
    plt.title('Iceberg: {}\nStart time: {}'.format(iip_berg.id_num, iip_berg.times[0]))
    i = 0
    filename = './drift_track_{}'.format(iip_berg.id_num)
    i = 0
    while True:
        i += 1
        newname = '{}_{:d}.png'.format(filename, i)
        if os.path.exists(newname):
            continue
        plt.savefig(newname)
        break
    plt.show()


def plot1(berg_df, iceberg):
    # t is the last timestep the model used
    t = iceberg[0,-1]
    t_all = iceberg[0,:]
    nbt = berg_df.loc[abs(berg_df['hrs_since'] - t) == min(abs(berg_df['hrs_since'] -t))].index.values[0] - berg_df.index[0] # nearest berg time index from 0
    plt.plot(berg_df['LONGITUDE'].values[:nbt+1], berg_df['LATITUDE'].values[:nbt+1], 
             label='observed', marker='o') 
    plt.plot(iceberg[1,:],iceberg[2,:], label='computed')
    for k in range(nbt+1):
        plt.annotate('{0:.1f}'.format(berg_df['hrs_since'].values[k]), 
                     xy=(berg_df['LONGITUDE'].values[k],
                         berg_df['LATITUDE'].values[k]))
        tmp_diff = abs(t_all - berg_df['hrs_since'].values[k])
        tmp_i = np.where(tmp_diff == tmp_diff.min())[0][0]
        if tmp_i >= iceberg[0,:].size:
            plt.annotate('{0:.1f}'.format(iceberg[0,-1]), 
                     xy=(iceberg[1,-1],iceberg[2,-1]))
            plt.plot(iceberg[1,-1],iceberg[2,-1], marker='o', color='black')
        else:
            plt.annotate('{0:.1f}'.format(t_all[tmp_i]), 
                     xy=(iceberg[1,tmp_i],iceberg[2,tmp_i]))
            plt.plot(iceberg[1,tmp_i],iceberg[2,tmp_i], marker='o', color='black')

        
    plt.plot(iceberg[1,:],iceberg[2,:], label='computed')
    plt.xlabel('Longitude'); plt.ylabel('Latitude')
    plt.title('Iceberg Drift\nSeason year: {}\nIIP indices: {}'.format(season_year, chosen_inds_arr[chosen_track_ind][:nbt+1]))
    plt.savefig('./drift_track1.png')
    plt.show()


def plot2(berg_df, iceberg):
    
    obs_lats = berg_df['LATITUDE'].tolist()
    obs_lons = berg_df['LONGITUDE'].tolist()
    obs_times = berg_df['hrs_since'].tolist()
    com_lats = iceberg[2,:]
    com_lons = iceberg[1,:]
    com_times = iceberg[0,:]
    
    # find the max time (in hours)
    max_t = min(max(obs_times), max(com_times))
    
    # trim data to max time
    obs_ind = bisect.bisect_right(obs_times, max_t)
    obs_times = obs_times[:obs_ind+1]
    obs_lats = obs_lats[:obs_ind+1]
    obs_lons = obs_lons[:obs_ind+1]
    
    # find a bounding box in space
    min_lat = min(min(obs_lats), min(com_lats))
    max_lat = max(max(obs_lats), max(com_lats))
    min_lon = min(min(obs_lons), min(com_lons))
    max_lon = max(max(obs_lons), max(com_lons))
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    buf = 0.1
    ax.set_extent([min_lon-buf, max_lon+buf, min_lat-buf, max_lat+buf], ccrs.PlateCarree())
    #ax.stock_img()
    ax.coastlines('50m')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    
    
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
        
    labels = [round(elem, 1) for elem in obs_times] 
    for label, x, y in zip(labels, obs_lons, obs_lats):
        ax.plot(x, y, marker='o', color='red')
        ax.text(x, y, label)
    
    for i, i2 in enumerate(com_times):
        for j in labels:
            if abs(com_times[i] - j) < 0.001:
                ax.plot(com_lons[i], com_lats[i], marker='X', color='green')
                ax.text(com_lons[i], com_lats[i], j)
    
    ax.plot(obs_lons, obs_lats, color='red')
    ax.plot(com_lons, com_lats, color='green')
    
    plt.savefig('./drift_track2.png')
    plt.show()


def plot3(berg_df, iceberg):

    obs_lats = berg_df['LATITUDE'].tolist()
    obs_lons = berg_df['LONGITUDE'].tolist()
    obs_times = berg_df['hrs_since'].tolist()
    com_lats = iceberg[2,:]
    com_lons = iceberg[1,:]
    com_times = iceberg[0,:]
    
    # find the max time (in hours)
    max_t = min(max(obs_times), max(com_times))
    
    # trim data to max time
    obs_ind = bisect.bisect_right(obs_times, max_t)
    obs_times = obs_times[:obs_ind+1]
    obs_lats = obs_lats[:obs_ind+1]
    obs_lons = obs_lons[:obs_ind+1]
    
    # find a bounding box in space
    min_lat = min(min(obs_lats), min(com_lats))
    max_lat = max(max(obs_lats), max(com_lats))
    min_lon = min(min(obs_lons), min(com_lons))
    max_lon = max(max(obs_lons), max(com_lons))
    
    m = Basemap(llcrnrlon=min_lon-0.1, llcrnrlat=min_lat-0.1, urcrnrlon=max_lon+0.1, urcrnrlat=max_lat+0.1)
    
    m.drawmapboundary(fill_color='aqua')
    m.fillcontinents(color='coral',lake_color='aqua')

    # labels = [left,right,top,bottom]
    parallels = np.arange(min_lat, max_lat, 0.1)
    m.drawparallels(parallels,labels=[False,True,False,False])
    meridians = np.arange(min_lon, max_lon, 0.1)
    meridians = m.drawmeridians(meridians,labels=[False,False,True,False])
    for mer in meridians:
        try:
            meridians[mer][1][0].set_rotation(45)
        except:
            pass
    
    obs_xpt, obs_ypt = m(obs_lons, obs_lats)
    com_xpt, com_ypt = m(com_lons, com_lats)
    m.plot(obs_xpt, obs_ypt, 'o')
    m.plot(com_xpt, com_ypt, '.') 
    
    labels = [round(elem, 1) for elem in obs_times] 
    for label, x, y in zip(labels, obs_xpt, obs_ypt):
        plt.text(x, y, label)
    
    for i, i2 in enumerate(com_times):
        for j in labels:
            if abs(com_times[i] - j) < 0.001:
                plt.plot(com_lons[i], com_lats[i], marker='X', color='black')
                plt.text(com_lons[i], com_lats[i], j)
    
    #m.drawmapscale(lonn-buf/1.5, lat0+buf/5, lon0, lat0, length=100, barstyle='fancy')
    plt.xlabel('Latitude')
    plt.ylabel('Longitude')
    plt.savefig('./drift_track3.png')
    plt.show()



def plot_return(iip_berg, mod_berg):
    f = plt.figure()
    tol = 0.1  # tstep
    for i,t in enumerate(iip_berg.t2000):
        diff = abs(np.asarray(mod_berg.t2000) - t)
        for j,diff_t in enumerate(diff):
            if diff_t < tol:
                break
        plt.plot(iip_berg.lons[i], iip_berg.lats[i], marker='o', color='red')
        plt.text(iip_berg.lons[i], iip_berg.lats[i], '{0:.1f}'.format(t-iip_berg.t2000[0]))
        plt.plot(mod_berg.lons[j], mod_berg.lats[j], marker='o', color='yellow')
        plt.text(mod_berg.lons[j], mod_berg.lats[j], '{0:.1f}'.format(t-iip_berg.t2000[0]))
    plt.plot(iip_berg.lons, iip_berg.lats, label='observed', color='red')
    plt.plot(mod_berg.lons, mod_berg.lats, label='computed', color='orange')
    plt.legend()
    plt.xlabel('Longitude'); plt.ylabel('Latitude')
    plt.title('Iceberg: {}\nStart time: {}'.format(iip_berg.id_num, iip_berg.times[0]))
    return f


def plot_return_size_vary(iip_berg, mod_berg_gr, mod_berg_bb, mod_berg_sm, mod_berg_med, 
                         mod_berg_lg, mod_berg_vlg, ind=None):
    f = plt.figure()
    tol = 0.1  # tstep
    for i,t in enumerate(iip_berg.t2000):
        diff = abs(np.asarray(mod_berg_med.t2000) - t)
        for j,diff_t in enumerate(diff):
            if diff_t < tol:
                break
        plt.plot(iip_berg.lons[i], iip_berg.lats[i], marker='o', color='red')
        plt.text(iip_berg.lons[i], iip_berg.lats[i], '{0:.1f}'.format(t-iip_berg.t2000[0]))
        plt.plot(mod_berg_med.lons[j], mod_berg_med.lats[j], marker='o', color='yellow')
        plt.text(mod_berg_med.lons[j], mod_berg_med.lats[j], '{0:.1f}'.format(t-iip_berg.t2000[0]))
    plt.plot(iip_berg.lons, iip_berg.lats, label='observed', color='red')
    plt.plot(mod_berg_gr.lons, mod_berg_gr.lats, label='gr', color='orange')
    plt.plot(mod_berg_bb.lons, mod_berg_bb.lats, label='bb', color='green')
    plt.plot(mod_berg_sm.lons, mod_berg_sm.lats, label='sm', color='blue')
    plt.plot(mod_berg_med.lons, mod_berg_med.lats, label='med', color='black')
    plt.plot(mod_berg_lg.lons, mod_berg_lg.lats, label='lg', color='purple')
    plt.plot(mod_berg_vlg.lons, mod_berg_vlg.lats, label='vlg', color='yellow')
    plt.legend()
    plt.xlabel('Longitude'); plt.ylabel('Latitude')
    if ind is not None:
        plt.title('Index: {}, Iceberg: {}\nStart time: {}'.format(ind,iip_berg.id_num, iip_berg.times[0]))
    else:
        plt.title('Iceberg: {}\nStart time: {}'.format(iip_berg.id_num, iip_berg.times[0]))
    return f


def plot_return_size_vary_no_time(iip_berg, mod_berg_gr, mod_berg_bb, mod_berg_sm, mod_berg_med, 
                         mod_berg_lg, mod_berg_vlg, ind=None):
    f = plt.figure()

    plt.plot(iip_berg.lons, iip_berg.lats, label='observed', color='red')
    plt.plot(mod_berg_gr.lons, mod_berg_gr.lats, label='gr', color='orange')
    plt.plot(mod_berg_bb.lons, mod_berg_bb.lats, label='bb', color='green')
    plt.plot(mod_berg_sm.lons, mod_berg_sm.lats, label='sm', color='blue')
    plt.plot(mod_berg_med.lons, mod_berg_med.lats, label='med', color='black')
    plt.plot(mod_berg_lg.lons, mod_berg_lg.lats, label='lg', color='purple')
    plt.plot(mod_berg_vlg.lons, mod_berg_vlg.lats, label='vlg', color='yellow')
    plt.legend()
    plt.xlabel('Longitude'); plt.ylabel('Latitude')

    return f


def animate_winds(atm_data, iip_berg, mod_berg):

    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([atm_data.x_min, atm_data.x_max, atm_data.y_min, atm_data.y_max], ccrs.PlateCarree())
    #ax.stock_img()
    ax.coastlines('50m')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ax.plot(iip_berg.lons[:], iip_berg.lats[:], color='orange')
    ax.plot(mod_berg.lons[:], mod_berg.lats[:], color='yellow')
    wind_mag = np.sqrt(atm_data.UA**2 + atm_data.VA**2)
    im = plt.imshow(wind_mag[0,:,:], 
                    extent=[atm_data.lons[0], atm_data.lons[-1] + atm_data.xy_res, atm_data.lats[0], atm_data.lats[-1] + atm_data.xy_res],
                    origin = 'lower', vmin=2, vmax=8)
    plt.colorbar()
    quiv = plt.quiver(atm_data.lons[:] + atm_data.xy_res, atm_data.lats[:] + atm_data.xy_res, atm_data.UA[0,:,:].T, atm_data.VA[0,:,:].T, 
                      scale=20, headwidth=5, width=0.005)
    title = plt.title('')


    def animate(i):

        im.set_data(wind_mag[i,:,:])
        quiv.set_UVC(atm_data.UA[i,:,:].T, atm_data.VA[i,:,:].T)
        title.set_text('time: {:.0f} hours'.format(i*6))
        return im

    anim = FuncAnimation(fig, animate, frames=wind_mag[:,0,0].size-1, interval=1000)
    #HTML(anim.to_html5_video())
    anim.save('plots/wind_mag.gif',writer='imagemagick')


def animate_currents(ocean_data, iip_berg, mod_berg):

    fig = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([ocean_data.x_min, ocean_data.x_max, ocean_data.y_min, ocean_data.y_max], ccrs.PlateCarree())
    #ax.stock_img()
    ax.coastlines('50m')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ax.plot(iip_berg.lons[:], iip_berg.lats[:], color='yellow')
    ax.plot(mod_berg.lons[:], mod_berg.lats[:], color='red')
    water_mag = np.sqrt(ocean_data.UW**2 + ocean_data.VW**2)

    im = plt.imshow(water_mag[0,:,:], 
                    extent=[ocean_data.lons[0], ocean_data.lons[-1] + ocean_data.xy_res, ocean_data.lats[0], ocean_data.lats[-1] + ocean_data.xy_res],
                    origin = 'lower',vmin=0, vmax=0.3)
    plt.colorbar()
    quiv = plt.quiver(ocean_data.lons[:] + ocean_data.xy_res, ocean_data.lats[:] + ocean_data.xy_res, ocean_data.UW[0,:,:].T, ocean_data.VW[0,:,:].T, 
                      scale=1, headwidth=5, width=0.005)
    title = plt.title('')

    def animate(i):

        im.set_data(water_mag[i,:,:])
        quiv.set_UVC(ocean_data.UW[i,:,:].T, ocean_data.VW[i,:,:].T)
        title.set_text('time: {:.0f} hours'.format(i))
        return im

    anim = FuncAnimation(fig, animate, frames=water_mag[:,0,0].size-1, interval=1000)
    #HTML(anim.to_html5_video())
    anim.save('plots/water_mag_9.gif',writer='imagemagick')
