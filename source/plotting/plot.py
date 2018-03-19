import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.basemap import Basemap
import bisect
from matplotlib.animation import FuncAnimation
import netCDF4 as nc


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



        
def animate_currents(ocean_data, iip_berg, mod_berg):
    
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
        quiv.set_UVC(UW[i,:,:].T, VW[i,:,:].T)
        title.set_text('time: {:.0f} hours'.format(mod_berg_t1950[i]-mod_berg_t1950[0]))
        line.set_data(mod_berg.lons[0:i+1], mod_berg.lats[0:i+1])
        
        
        return im, line
    
    anim = FuncAnimation(fig, animate, frames=water_mag[:,0,0].size-1, interval=100)
    #HTML(anim.to_html5_video())
    anim.save('plots/water_mag.gif',writer='imagemagick')
        
    
    
    
    
def animate_winds(atm_data, iip_berg, mod_berg):
    
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
    quiv = plt.quiver(atm_data.lons[lon0-1:lonn+1], atm_data.lats[lat0-1:latn+1], UA[0,:,:].T, VA[0,:,:].T, 
                      scale=20, headwidth=5, width=0.005)
    title = plt.title('time: 0 hours')


    def animate(i):

        im.set_data(wind_mag[i,:,:])
        quiv.set_UVC(UA[i,:,:].T, VA[i,:,:].T)
        title.set_text('time: {:.0f} hours'.format(mod_berg_t1900[i]-mod_berg_t1900[0]))
        line.set_data(mod_berg.lons[0:i+1], mod_berg.lats[0:i+1])
        
        
        return im, line
    
    anim = FuncAnimation(fig, animate, frames=wind_mag[:,0,0].size-1, interval=100)
    #HTML(anim.to_html5_video())
    anim.save('plots/wind_mag.gif',writer='imagemagick')
    

def plot_turnbull(iip_berg, mod_berg):
    
    f = plt.figure()
    
    iip_times = np.asarray(iip_berg.datetimes)
    mod_times = np.asarray(mod_berg.datetimes)
    iip_t0 = iip_berg.datetimes[0]
    mod_t0 = mod_berg.datetimes[0]
    assert iip_t0 == mod_t0
    iip_times0 = iip_times - iip_t0
    mod_times0 = mod_times - mod_t0
    iip_hours0, mod_hours0 = np.empty(0), np.empty(0)
    
    for i in range(len(iip_times0)):
        iip_hours0 = np.append(iip_hours0, round(iip_times0[i].days*24 + iip_times0[i].seconds/3600, 1))
        plt.text(iip_berg.lons[i], iip_berg.lats[i], '{}'.format(iip_hours0[i]))
        
    for j in range(len(mod_times0)):
        mod_hours0 = np.append(mod_hours0,round(mod_times0[j].days*24 + mod_times0[j].seconds/3600, 1))
        #plt.text(mod_berg.lons[j], mod_berg.lats[j], '{}H'.format(mod_hours0[j]))
    ind_arr = []
    tol = 1e-3
    for k in range(len(iip_times0)):
        mod_diff = mod_hours0 - iip_hours0[k]
        #ind = np.where(mod_diff == 0.0)[0][0]
        ind = np.where(abs(mod_diff < tol))[0][0]
        ind_arr.append(ind)
        plt.scatter(mod_berg.lons[ind], mod_berg.lats[ind], color='orange')
        plt.text(mod_berg.lons[ind], mod_berg.lats[ind], '{}'.format(mod_hours0[ind]))
        
    
        
    plt.scatter(iip_berg.lons, iip_berg.lats, label='observed', color='red')
    plt.plot(mod_berg.lons, mod_berg.lats, label='computed', color='orange')

    plt.legend()
    plt.xlabel('Longitude'); plt.ylabel('Latitude')

    return f


def plot_turnbull_subplots(iip_berg, mod_berg):
    
    f, (ax1, ax2) = plt.subplots(1,2)
    f.tight_layout()
    
    iip_times = np.asarray(iip_berg.datetimes)
    mod_times = np.asarray(mod_berg.datetimes)
    iip_t0 = iip_berg.datetimes[0]
    mod_t0 = mod_berg.datetimes[0]
    assert iip_t0 == mod_t0
    iip_times0 = iip_times - iip_t0
    mod_times0 = mod_times - mod_t0
    iip_hours0, mod_hours0 = np.empty(0), np.empty(0)
    
    for i in range(len(iip_times0)):
        iip_hours0 = np.append(iip_hours0, round(iip_times0[i].days*24 + iip_times0[i].seconds/3600, 1))
        ax1.text(iip_berg.lons[i], iip_berg.lats[i], '{}'.format(iip_hours0[i]))
        
    for j in range(len(mod_times0)):
        mod_hours0 = np.append(mod_hours0,round(mod_times0[j].days*24 + mod_times0[j].seconds/3600, 1))
        #plt.text(mod_berg.lons[j], mod_berg.lats[j], '{}H'.format(mod_hours0[j]))
    ind_arr = []
    for k in range(len(iip_times0)):
        mod_diff = mod_hours0 - iip_hours0[k]
        ind = np.where(mod_diff == 0.0)[0][0]
        ind_arr.append(ind)
        ax1.scatter(mod_berg.lons[ind], mod_berg.lats[ind], color='orange')
        ax1.text(mod_berg.lons[ind], mod_berg.lats[ind], '{}'.format(mod_hours0[ind]))
        
    ax1.set_xlabel('Longitude') 
    ax2.set_ylabel('Latitude')  
    ax1.scatter(iip_berg.lons, iip_berg.lats, label='observed', color='red')
    ax1.plot(mod_berg.lons, mod_berg.lats, label='computed', color='orange')

    dists = []
    dist_times = []
    for l, item in enumerate(ind_arr):
        dist = np.sqrt( (mod_berg.lons[item] - iip_berg.lons[l])**2 + \
                       (mod_berg.lats[item] - iip_berg.lats[l])**2 )
        dists.append(dist)
        dist_time = mod_hours0[item]
        dist_times.append(dist_time)
    
    ax2.set_ylabel('Distance (deg)')
    ax2.set_xlabel('Time (hrs)')
    ax2.set_ylim(0, max(dists) + np.mean(dists))
    ax2.scatter(dist_times, dists)

    return f