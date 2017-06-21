from iceberg import Iceberg
from config import *
from functions import *
from store_objects import *
from load_objects import *
from analysis import *
from test import *
import numpy as np
from scipy.interpolate import interpn
from matplotlib.backends.backend_pdf import PdfPages


def main():
    global bb
    plot_list = []
    save_plots = True
    for bb in bvec:
        print("Iceberg size class: {}".format(bb))
        silent_remove('bergClass{}.pkl'.format(bb))
        L, W, H = bergdims[bb-1,0],bergdims[bb-1,1],bergdims[bb-1,2]
        VOL = L*W*H; dL,dW,dH,dVOL = 0,0,0,0
        global j
        for j in range(0,trajnum):
            assert_tol_matrix(LAT,np.ravel(mLAT),0,j)
            assert_tol_matrix(LON,np.ravel(mLON),0,j)
            assert_tol_matrix(msk,mmsk,0,j)
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
        plot_name = 'plot' + str(bb)
        plot_name = plot_model_diff(dXIL,dYIL)
        plot_list.append(plot_name)
    
    if save_plots:
        with PdfPages('plots.pdf') as pdf:
            for plot in plot_list:
                pdf.savefig(plot)

        
def drift(I,loc,dims):

    OB = False
    interpolate = False

    timestep = tt[tts + I]                                      
    t=timestep
    t1  = int(np.floor(timestep)); t2 = t1 + 1 

    if interpolate:

        ti1=(t2-t)/(t2-t1)
        ti2=(t-t1)/(t2-t1)
        ti=t1*ti1+t2*ti2

        XI1 = np.where(LON <= loc[I,0])[0][-1]                                    
        XI2 = np.where(LON > loc[I,0])[0][0]                                      
        YI1 = np.where(LAT <= loc[I,1])[0][-1]                                    
        YI2 = np.where(LAT > loc[I,1])[0][0] 

        points = ((XI1,XI2),(YI1,YI2),(t1,t2))

        x=loc[I,0]; 
        xi1=(LON[XI2]-x)/(LON[XI2]-LON[XI1]); xi2=(x-LON[XI1])/(LON[XI2]-LON[XI1])
        xi=XI1*xi1+XI2*xi2
        
        y=loc[I,1]
        yi1=(LAT[YI2]-y)/(LAT[YI2]-LAT[YI1]); yi2=(y-LAT[YI1])/(LAT[YI2]-LAT[YI1])
        yi=YI1*yi1+YI2*yi2
        
        xyti = [xi,yi,ti]

        fields = [uaF,vaF,uwF,vwF,sst]; name = []
        ii = 0
        for field in fields:
            values = field[XI1:XI2+1,YI1:YI2+1,t1:t2+1]
            name.append(interpn(points,values,xyti)[0]) 
            ii+=1
        ua = name[0]; va = name[1]; uw = name[2]; vw = name[3]; SST = name[4]

    else:

        # Find nearest neighbour
        XI = find_nearest(LON, loc[I,0])
        YI = find_nearest(LAT, loc[I,1])
        assert_tol(XI,mXI[bb-1,j,I]-1,I,j)     
        assert_tol(YI,mYI[bb-1,j,I]-1,I,j)     

        # Interpolate fields linearly between timesteps
        dt1 = timestep - t1; dt2 = t2 - timestep
        
        ua = uaF[XI,YI,t1] * dt1 + uaF[XI,YI,t2] * dt2 
        va = vaF[XI,YI,t1] * dt1 + vaF[XI,YI,t2] * dt2 
        uw = uwF[XI,YI,t1] * dt1 + uwF[XI,YI,t2] * dt2 
        vw = vwF[XI,YI,t1] * dt1 + vwF[XI,YI,t2] * dt2 
        SST = sst[XI,YI,t1] * dt1 + sst[XI,YI,t2] * dt2
        assert_tol(ua,mUA[bb-1,j,I],I,j)
        assert_tol(va,mVA[bb-1,j,I],I,j)
        assert_tol(uw,mUW[bb-1,j,I],I,j)
        assert_tol(vw,mVW[bb-1,j,I],I,j)

    # Compute wind speed and "U tilde" at location for a given icesize
    Ua = np.sqrt(ua**2 + va**2)
    UT = Ut(Ua,loc[I,1], S(dims[I,0],dims[I,1]),Cw,g,om)

    # now compute analytic icevelocity solution
    ui = uw-g*a(UT)*va+g*b(UT)*ua
    vi = vw+g*a(UT)*ua+g*b(UT)*va
    assert_tol(ui,mUI[bb-1,j,I],I,j)
    assert_tol(vi,mVI[bb-1,j,I],I,j)

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

    assert_tol(loc[I+1,0],mXIL[bb-1,j,I+1],I,j)
    assert_tol(loc[I+1,1],mYIL[bb-1,j,I+1],I,j)

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

    else:
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

    #assert_tol(dims[I+1,0],mL[bb-1,j,I+1],I,j,correct=True,tol=1e-9)
    #assert_tol(dims[I+1,1],mW[bb-1,j,I+1],I,j,correct=True,tol=1e-9)
    #assert_tol(dims[I+1,2],mH[bb-1,j,I+1],I,j,correct=True,tol=1e-9)
    #assert_tol(dims[I+1,3],mVOL[bb-1,j,I+1],I,j,correct=True)
    #assert_tol(ddims[I+1,3],mDVOL[bb-1,j,I+1],I,j,correct=True)

    return dims,ddims,melted


if __name__=="__main__":
    main()
