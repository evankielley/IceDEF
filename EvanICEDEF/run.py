from iceberg import Iceberg
from constants import *
from functions import *
from store_objects import *
from load_objects import *
from analysis import *

def main():

    for bb in bvec:
        silent_remove('bergClass{}.pkl'.format(bb))
        L, W, H = bergdims[bb-1,0]*1.0,bergdims[bb-1,1]*1.0,bergdims[bb-1,2]*1.0
        VOL = L*W*H; dL,dW,dH,dVOL = 0,0,0,0
        for j in range(0,trajnum):
            berg = Iceberg(nt)
            berg.dims[0,:] = L,W,H,VOL
            berg.dimsChange[0,:] = dL,dW,dH,dVOL
            #ts = np.random.randint(0,round(startrange)) 
            ts = mts[bb-1,j]
            global tts
            tts = ts*tres
            lt = nt-tts
            #xig = seed_X[np.random.randint(1, seed_X.size)] 
            #yig = seed_Y[np.random.randint(1, seed_Y.size)]
            randoX = mrandoX[bb-1,j]; randoY = mrandoY[bb-1,j]
            xig = seed_X[randoX-1]; yig = seed_Y[randoY-1] 
            berg.location[0,:] = LON[xig-1], LAT[yig-1]
            i=0
            while not berg.outOfBounds and not berg.melted and i<lt-1:
                I=i
                berg.location,berg.outOfBounds,Ua,SST,ui,uw,vi,vw = drift(I,berg.location,berg.dims)
                berg.dims,berg.dimsChange,berg.melted = melt(I,berg.dims,berg.dimsChange,Ua,SST,ui,uw,vi,vw)
                i += 1
            store_objects(berg, 'bergClass{}.pkl'.format(bb))
        XIL,YIL = load_objects(pyOutloc + 'bergClass{}.pkl'.format(bb),trajnum,nt)
        dXIL, dYIL = compare_outputs(mXIL,mYIL,bb,pyOutloc+'bergClass{}.pkl'.format(bb),trajnum,nt)
        plot_model_diff(dXIL,dYIL)


def drift(I,loc,dims):

    OB = False

    # Find nearest neighbour
    XI = int(find_nearest(LON, loc[I,0]))
    YI = int(find_nearest(LAT, loc[I,1]))
    
    # Interpolate fields linearly between timesteps
    timestep = tt[tts + I]
    t1  = int(np.floor(timestep)); t2 = t1 + 1   
    dt1 = timestep - t1; dt2 = t2 - timestep
    
    ua = uaF[XI,YI,t1] * dt1 + uaF[XI,YI,t2] * dt2 
    va = vaF[XI,YI,t1] * dt1 + vaF[XI,YI,t2] * dt2 
    uw = uwF[XI,YI,t1] * dt1 + uwF[XI,YI,t2] * dt2 
    vw = vwF[XI,YI,t1] * dt1 + vwF[XI,YI,t2] * dt2 
    SST = sst[XI,YI,t1] * dt1 + sst[XI,YI,t2] * dt2

    # Compute wind speed and "U tilde" at location for a given icesize
    Ua = np.sqrt(ua**2 + va**2)
    UT = Ut(Ua,loc[I,1], S(dims[I,0],dims[I,1]),Cw,g,om)

    # now compute analytic icevelocity solution
    ui = uw-g*a(UT)*va+g*b(UT)*ua
    vi = vw+g*a(UT)*ua+g*b(UT)*va

    # Icetranslation -- Note the conversion from meters to degrees lon/lat   
    dlon = ui*dtR 
    dlat = vi*dtR
    loc[I+1,1] = loc[I,1] + dlat
    loc[I+1,0] = loc[I,0]+ dlon/np.cos((loc[I+1,1]+loc[I,1])/2*np.pi/180)

    # Check if out-of-bounds
    if loc[I+1,0]>max(LON) or loc[I+1,0]<min(LON) or loc[I+1,1]>max(LAT) or loc[I+1,1]<min(LAT):
        OB = True
    else:
        xi2a,xi2b,yi2a,yi2b = find_berg_grid(LAT,LON,loc[I+1,0], loc[I+1,1])
        if msk[xi2a,yi2a]==0 or msk[xi2a,yi2b]==0 or msk[xi2b,yi2a]==0 or msk[xi2b,yi2b]==0:
            loc[I+1,0] = loc[I,0]
            loc[I+1,1] = loc[I,1]

    return loc,OB,Ua,SST,ui,uw,vi,vw


def melt(I,dims,ddims,Ua,SST,ui,uw,vi,vw):
    melted = False
    # Melt terms
    Me = CMe1 * (Cs1 * Ua**Cs2 + Cs3 * Ua)  ## Wind driven erosion
    Mv = CMv1 * SST + CMv2 * SST**2  ## Thermal side wall erosion 
    Mb = CMb1*np.power(np.sqrt(np.square(ui-uw)+np.square(vi-vw)),CMb2)*(SST-Ti0)/dims[I,0]**CMb3

    # Melt rates
    ddims[I,0] = - Mv - Me 
    ddims[I,1] = - Mv - Me 
    ddims[I,2] = - Mb
    dims[I+1,0] = dims[I,0]+ddims[I,0]*Dt 
    dims[I+1,1] = dims[I,1]+ddims[I,1]*Dt 
    dims[I+1,2] = dims[I,2]+ddims[I,2]*Dt 

    # Check if iceberg size is negative
    if dims[I+1,0]<0 or dims[I+1,1]<0 or dims[I+1,2]<0:
        dims[I+1,0] = 0; dims[I+1,1] = 0; dims[I+1,2] = 0
        melted = True

    # Rollover
    if dims[I+1,1] < (0.85*dims[I+1,2]):
        hn = dims[I+1,1]  ## new height
        dims[I+1,1] = dims[I+1,2] 
        dims[I+1,2] = hn

    # Check if length is greater than width
    if dims[I+1,1] > dims[I+1,0]:
        wn = dims[I+1,0] 
        dims[I+1,0] = dims[I+1,1] 
        dims[I+1,1] = wn

    # New volume and change in volume (dv)
    dims[I+1,3] = dims[I+1,0]*dims[I+1,1]*dims[I+1,2]
    ddims[I+1,3] = dims[I,3] - dims[I+1,3]

    # Store melt rates
    Mev = Me; Mvv = Mv; Mbv = Mb            

    return dims,ddims,melted


if __name__=="__main__":
    main()
