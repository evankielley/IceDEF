from iceberg import Iceberg
from load_inputs import *
from param_funcs import *
from find_neighbours import *
from make_plots import *
from store_objects import *
from load_objects import *
from compare_outputs import *
from test import *
from load_mat_vars import *
from analysis import *

# INPUTS ##############################################################
# Paths
modelfull,modelshort,root,condloc,outloc,modelloc,pyOutloc = load_paths()
# Constants
R,rhow,rhoa,rhoi,drho,Cw,Ca,om,g = load_constants()
# Melt parameters
Ti0,Cs1,Cs2,Cs3,CMv1,CMv2,CMe1,CMb1,CMb2,CMb3 = load_melt_params()
# Run parameters
bvec,trajnum,final_t,startrange,tres,DT,Dt,dt,R,dtR,t,nt,tt = load_run_params()
# Model inputs
msk,bergdims,seed_X,seed_Y,LAT,LON,uwF,vwF,uaF,vaF,sst = load_fields(condloc,modelshort,modelloc,t)
# Fixed inputs
fixed_file_path = outloc + 'fixed.mat'
mts, mrandoX, mrandoY = load_fixed_vars(fixed_file_path)
# Expected values
expected_file_path = outloc + 'output_full.mat'
mXIL,mYIL,mVOL,mDVOL,mUI,mVI,mUA,mVA,mUW,mVW,mTE,mMemat,mMvmat,mMbmat = load_expected_vars(expected_file_path)
######################################################################

corrector = False 
ignorer = True 

for bb in bvec:

    silent_remove('bergClass{}.pkl'.format(bb))

    L, W, H = bergdims[bb-1,0]*1.0,bergdims[bb-1,1]*1.0,bergdims[bb-1,2]*1.0
    VOL = L*W*H
    dL,dW,dH,dVOL = 0,0,0,0
    
    for j in range(0,trajnum):
        
        berg = Iceberg(nt)
        berg.dims[0,:] = L,W,H,VOL
        berg.dimsChange[0,:] = dL,dW,dH,dVOL
        
        #ts = np.random.randint(0,round(startrange)) 
        ts = mts[bb-1,j]
        tts = ts*tres
        lt = nt-tts

        #xig = seed_X[np.random.randint(1, seed_X.size)] 
        #yig = seed_Y[np.random.randint(1, seed_Y.size)]
        randoX = mrandoX[bb-1,j]; randoY = mrandoY[bb-1,j]
        xig = seed_X[randoX-1]; yig = seed_Y[randoY-1] 

        berg.location[0,:] = LON[xig-1], LAT[yig-1]

        i=0;z=0
        while not berg.outOfBounds and not berg.melted and i<lt-1:
            I=i
# Drift
            # Find nearest neighbour
            XI = find_nearest(LON, berg.location[I,0])+1
            YI = find_nearest(LAT, berg.location[I,1])+1
            #XI = mXI[z]
            #YI = mYI[z]
            
            #assert_tol(XI,mXI[z],z)
            #assert_tol(YI,mYI[z],z)
            
            # Interpolate fields linearly between timesteps
            timestep = tt[tts + I]+1
            #timestep = mtimestep[z]
            
            t1  = np.floor(timestep); t2 = t1 + 1   
            #t1 = mt1[z]; t2 = t1+1

            #assert_tol(t1,mt1[z],z)    
    
            dt1 = timestep - t1; dt2 = t2 - timestep
            XI = int(XI)-1; YI = int(YI)-1 
            t1 = int(t1)-1; t2 = int(t2)-1
            
            ua = uaF[XI,YI,t1] * dt1 + uaF[XI,YI,t2] * dt2 
            va = vaF[XI,YI,t1] * dt1 + vaF[XI,YI,t2] * dt2 
            uw = uwF[XI,YI,t1] * dt1 + uwF[XI,YI,t2] * dt2 
            vw = vwF[XI,YI,t1] * dt1 + vwF[XI,YI,t2] * dt2 
            SST = sst[XI,YI,t1] * dt1 + sst[XI,YI,t2] * dt2

            ua = assert_tol(ua,mUA[bb-1,j,z],z,corrector,ignorer)
            va = assert_tol(va,mVA[bb-1,j,z],z,corrector,ignorer)
            uw = assert_tol(uw,mUW[bb-1,j,z],z,corrector,ignorer)
            vw = assert_tol(vw,mVW[bb-1,j,z],z,corrector,ignorer)
            #assert_tol(SST,mSST[z],z)

            #ua = mua[z]; va = mva[z]; uw = muw[z]; vw = mvw[z]; SST = mSST[z]
            # Compute wind speed and "U tilde" at location for a given iceberg size
            Ua = np.sqrt(ua**2 + va**2)
            UT = Ut(Ua,berg.location[I,1], S(berg.dims[I,0],berg.dims[I,1]),Cw,g,om)

            #assert_tol(Ua,mUa[z],z)            
            berg.location[I,1] = assert_tol(berg.location[I,1],mYIL[bb-1,j,z],z,corrector,ignorer)
            berg.location[I,0] = assert_tol(berg.location[I,0],mXIL[bb-1,j,z],z,corrector,ignorer)
            #assert_tol(UT,mUT[z],z)

            # now compute analytic iceberg velocity solution
            ui = uw-g*a(UT)*va+g*b(UT)*ua
            vi = vw+g*a(UT)*ua+g*b(UT)*va

            ui = assert_tol(ui,mUI[bb-1,j,z],z,corrector,ignorer)
            vi = assert_tol(vi,mVI[bb-1,j,z],z,corrector,ignorer)

            #ui = mui[z]; vi = mvi[z]

            # Iceberg translation -- Note the conversion from meters to degrees lon/lat   
            dlon = ui * dtR 
            dlat = vi * dtR
            berg.location[I+1,1] = berg.location[I,1] + dlat
            berg.location[I+1,0] = berg.location[I,0]+ dlon/np.cos((berg.location[I+1,1]+berg.location[I,1])/2*np.pi/180)

            # Check if out-of-bounds
            if berg.location[I+1,0]>max(LON) or berg.location[I+1,0]<min(LON) or berg.location[I+1,1]>max(LAT) or berg.location[I+1,1]<min(LAT):
                berg.outOfBounds = True
            else:
                #xi2,yi2 = find_berg_grid(LAT,LON,berg.location[I+1,0], berg.location[I+1,1])
                xi2a,xi2b,yi2a,yi2b = find_berg_grid(LAT,LON,berg.location[I+1,0], berg.location[I+1,1])
                #if any(msk[xi2, yi2]) is 0:
                if msk[xi2a,yi2a]==0 or msk[xi2a,yi2b]==0 or msk[xi2b,yi2a]==0 or msk[xi2b,yi2b]==0:
                    berg.location[I+1,0] = berg.location[I,0]
                    berg.location[I+1,1] = berg.location[I,1]
            berg.location[I+1,1] = assert_tol(berg.location[I+1,1],mYIL[bb-1,j,z+1],z,corrector,ignorer)
            berg.location[I+1,0] = assert_tol(berg.location[I+1,0],mXIL[bb-1,j,z+1],z,corrector,ignorer)

# Melt
            # Melt terms
            Me = CMe1 * (Cs1 * Ua**Cs2 + Cs3 * Ua)  ## Wind driven erosion
            Mv = CMv1 * SST + CMv2 * SST**2  ## Thermal side wall erosion 
            Mb = CMb1*np.power(np.sqrt(np.square(ui-uw)+np.square(vi-vw)),CMb2)*(SST-Ti0)/berg.dims[I,0]**CMb3

            # Melt rates
            berg.dimsChange[I,0] = - Mv - Me 
            berg.dimsChange[I,1] = - Mv - Me 
            berg.dimsChange[I,2] = - Mb
            berg.dims[I+1,0] = berg.dims[I,0]+berg.dimsChange[I,0]*Dt 
            berg.dims[I+1,1] = berg.dims[I,1]+berg.dimsChange[I,1]*Dt 
            berg.dims[I+1,2] = berg.dims[I,2]+berg.dimsChange[I,2]*Dt 

            # Check if iceberg size is negative
            if berg.dims[I+1,0]<0 or berg.dims[I+1,1]<0 or berg.dims[I+1,2]<0:
                berg.dims[I+1,0] = 0; berg.dims[I+1,1] = 0; berg.dims[I+1,2] = 0
                berg.melted = True   ## Boolean

            # Rollover
            if berg.dims[I+1,1] < (0.85*berg.dims[I+1,2]):
                hn = berg.dims[I+1,1]  ## new height
                berg.dims[I+1,1] = berg.dims[I+1,2] 
                berg.dims[I+1,2] = hn

            # Check if length is greater than width
            if berg.dims[I+1,1] > berg.dims[I+1,0]:
                wn = berg.dims[I+1,0] 
                berg.dims[I+1,0] = berg.dims[I+1,1] 
                berg.dims[I+1,1] = wn

            # New volume and change in volume (dv)
            berg.dims[I+1,3] = berg.dims[I+1,0] * berg.dims[I+1,1] * berg.dims[I+1,2]
            berg.dimsChange[I+1,3] = berg.dims[I,3] - berg.dims[I+1,3]


            # Store melt rates
            Mev = Me; Mvv = Mv; Mbv = Mb            

            # Update iterator
            i += 1
            z += 1

        store_objects(berg, 'bergClass{}.pkl'.format(bb))

    XIL,YIL = load_objects(pyOutloc + 'bergClass{}.pkl'.format(bb),trajnum,nt)
    dXIL, dYIL = compare_outputs(mXIL,mYIL,bb,pyOutloc+'bergClass{}.pkl'.format(bb),trajnum,nt)
    make_plots(dXIL)
    make_plots(dYIL)
