import numpy as np
from iceberg import Iceberg
from specifications import *
from load_constants import *
from load_inputs import *
from load_objects import *
from functions import *
from store_objects import *
from analysis import *
from test import *
from scipy.interpolate import interpn
from matplotlib.backends.backend_pdf import PdfPages


def main():
    plot_list = []

    for bb in range(1, num_bergs+1):

        # Remove old output files
        silent_remove('berg{}.pkl'.format(bb))

        print("Iceberg: {}".format(bb))
        
        for j in range(0, num_trajs):

            # Create iceberg instance
            berg = Iceberg(nt, init_berg_dims[bb-1], init_berg_coords[bb-1])

            # Randomize start time
            ts = np.random.randint(0,round(final_t/2)) 
            tts = ts*tres
            lt = nt-tts

            # Simulate
            i=0
            while not berg.outOfBounds and not berg.melted and i<lt-1:
                I=i
                berg.coords,berg.outOfBounds,berg.grounded,Ua,SST,ui,uw,vi,vw = drift(I,tts,berg.coords,berg.dims)
                berg.dims,berg.melted = melt(I,berg.dims,Ua,SST,ui,uw,vi,vw)
                i += 1

            # Store iceberg instance
            store_objects(berg, 'berg{}.pkl'.format(bb))

        # Load all trajectories for one iceberg
        berg_lon, berg_lat = load_objects(root + 'berg{}.pkl'.format(bb), num_trajs, nt)
        
        # Plot
        plot_name = 'plot' + str(bb)
        plot_name = plot_track_on_map(berg_lon,berg_lat,show=True)
        plot_list.append(plot_name)
    
    # Save a PDF of all plots
    if save_plots:
        with PdfPages('plots.pdf') as pdf:
            for plot in plot_list:
                pdf.savefig(plot)

        
def drift(I,tts,berg_coords,berg_dims):

    GROUNDED = False
    OB = False

    timestep = tt[tts + I]                                      
    t = timestep
    t1 = int(np.floor(timestep)); t2 = t1 + 1 

    if interpolate:

        ti1=(t2-t)/(t2-t1)
        ti2=(t-t1)/(t2-t1)
        ti=t1*ti1+t2*ti2

        lon_1 = np.where(LON <= berg_coords[I,0])[0][-1]                                    
        lon_2 = np.where(LON > berg_coords[I,0])[0][0]                                      
        lat_1 = np.where(LAT <= berg_coords[I,1])[0][-1]                                    
        lat_2 = np.where(LAT > berg_coords[I,1])[0][0] 

        points = ((lon_1,lon_2),(lat_1,lat_2),(t1,t2))

        lon = berg_coords[I,0]; 
        lon1 = (LON[lon_2] - lon)/(LON[lon_2] - LON[lon_1]); 
        lon2 = (lon - LON[lon_1])/(LON[lon_2] - LON[lon_1])
        loni = lon_1*lon1 + lon_2*lon2
        
        lat = berg_coords[I,1]
        lat1 = (LAT[lat_2] - lat)/(LAT[lat_2] - LAT[lat_1]); 
        lat2 = (lat - LAT[lat_1])/(LAT[lat_2] - LAT[lat_1])
        lati = lat_1*lat1 + lat_2*lat2
        
        lon_lat_t_i = [loni,lati,ti]

        fields = [uaF,vaF,uwF,vwF,sst]; name = []
        ii = 0
        for field in fields:
            values = field[lon_1:lon_2+1,lat_1:lat_2+1,t1:t2+1]
            name.append(interpn(points,values,lon_lat_t_i)[0]) 
            ii+=1
        ua = name[0]; va = name[1]; uw = name[2]; vw = name[3]; SST = name[4]

    else:

        # Find nearest neighbour
        lon_ = find_nearest(LON, berg_coords[I,0])
        lat_ = find_nearest(LAT, berg_coords[I,1])

        # Interpolate fields linearly between timesteps
        dt1 = timestep - t1; dt2 = t2 - timestep
        
        ua = uaF[lon_,lat_,t1] * dt1 + uaF[lon_,lat_,t2] * dt2 
        va = vaF[lon_,lat_,t1] * dt1 + vaF[lon_,lat_,t2] * dt2 
        uw = uwF[lon_,lat_,t1] * dt1 + uwF[lon_,lat_,t2] * dt2 
        vw = vwF[lon_,lat_,t1] * dt1 + vwF[lon_,lat_,t2] * dt2 
        SST = sst[lon_,lat_,t1] * dt1 + sst[lon_,lat_,t2] * dt2

    # Compute wind speed and "U tilde" at location for a given icesize
    Ua = np.sqrt(ua**2 + va**2)
    UT = Ut(Ua,berg_coords[I,1], S(berg_dims[I,0],berg_dims[I,1]),Cw,g,om)

    # now compute analytic ice velocity solution
    ui = uw-g*alpha(UT)*va+g*beta(UT)*ua
    vi = vw+g*alpha(UT)*ua+g*beta(UT)*va

    # Icetranslation -- Note the conversion from meters to degrees lon/lat   
    dlon = ui*dtR 
    dlat = vi*dtR
    berg_coords[I+1,1] = berg_coords[I,1] + dlat
    berg_coords[I+1,0] = berg_coords[I,0]+ dlon/np.cos((berg_coords[I+1,1]+berg_coords[I,1])/2*np.pi/180)

    # Check if out-of-bounds
    if berg_coords[I+1,0]>max(LON) or berg_coords[I+1,0]<min(LON) or berg_coords[I+1,1]>max(LAT) or berg_coords[I+1,1]<min(LAT):
        OB = True
    else:
        lon2a,lon2b,lat2a,lat2b = find_berg_grid(LAT,LON,berg_coords[I+1,0], berg_coords[I+1,1])
        if msk[lon2a,lat2a]==0 or msk[lon2a,lat2b]==0 or msk[lon2b,lat2a]==0 or msk[lon2b,lat2b]==0:
            GROUNDED = True
            berg_coords[I+1,0] = berg_coords[I,0]
            berg_coords[I+1,1] = berg_coords[I,1]

    return berg_coords,OB,GROUNDED,Ua,SST,ui,uw,vi,vw


def melt(I,berg_dims,Ua,SST,ui,uw,vi,vw):

    MELTED = False

    # Melt terms
    Me = CMe1 * (Cs1 * Ua**Cs2 + Cs3 * Ua)  ## Wind driven erosion
    Mv = CMv1 * SST + CMv2 * SST**2  ## Thermal side wall erosion 
    Mb = CMb1*np.power(np.sqrt(np.square(ui-uw)+np.square(vi-vw)),CMb2)*(SST-Ti0)/berg_dims[I,0]**CMb3

    # New berg dims
    berg_dims[I+1,0] = berg_dims[I,0] - (Mv + Me)*Dt
    berg_dims[I+1,1] = berg_dims[I,1] - (Mv + Me)*Dt
    berg_dims[I+1,2] = berg_dims[I,2] - Mb*Dt

    # Check if iceberg size is negative
    if berg_dims[I+1,0]<0 or berg_dims[I+1,1]<0 or berg_dims[I+1,2]<0:
        berg_dims[I+1,0]=0; berg_dims[I+1,1]=0; berg_dims[I+1,2]=0
        MELTED = True

    else:
        # Rollover
        if berg_dims[I+1,1] < (0.85*berg_dims[I+1,2]):
            hn = berg_dims[I+1,1]  # new height
            berg_dims[I+1,1] = berg_dims[I+1,2] 
            berg_dims[I+1,2] = hn

        # Check if length is greater than width
        if berg_dims[I+1,1] > berg_dims[I+1,0]:
            wn = berg_dims[I+1,0] 
            berg_dims[I+1,0] = berg_dims[I+1,1] 
            berg_dims[I+1,1] = wn

    return berg_dims, MELTED


if __name__=="__main__":
    main()
