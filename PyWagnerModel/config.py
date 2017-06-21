import scipy.io as sio
import numpy as np
import numpy.matlib

R = np.float64(6378*1e3)  ## earth radius in m
rhow = 1027  ## density of water (kg/m^3)
rhoa = 1.2  ## density of air   (kg/m^3)
rhoi = 850  ## density of ice   (kg/m^3)
drho = rhow - rhoi
Cw = 0.9  ## bulk coefficient water  (Bigg et al 1997)
Ca = 1.3  ## bulk coefficient air    (Bigg et al 1997)
om = np.float64(7.2921e-5)  ## rotation rate of earth (rad/s)
g = np.sqrt(rhoa*drho/rhow/rhoi*(Ca/Cw))  ## gamma = np.sqrt(ca/cw)

Ti0 = -4
Cs1 = 1.5; Cs2 = 0.5; Cs3 = 0.1
CMv1 = 7.62e-3; CMv2 = 1.29e-3; CMe1 = 0.5
CMb1 = 0.58; CMb2 = 0.8; CMb3 = 0.2

bvec = range(1,11)
trajnum = 25            # total number of iceberg trajectories to compute
final_t = 122           # number of input field time steps
startrange = final_t/2  # input field start range
tres = 3                # time resoln such that "model Dt"="input DT"/tres
DT = 3                  # Input fields time step
Dt = DT/tres            # model timestep in days
dt = np.float64(Dt*24*3600)         # model timestep in seconds
R = np.float64(6378*1e3)            # earth radius in m
dtR = np.float64(dt/R*180/np.pi)       # need this ratio for distances in "drifting.m"
t = range(0, final_t)      # how long is the run
nt = len(t)*tres            # number of model timesteps
tt = np.linspace(0, len(t)-1,nt)  # model time

modelfull = 'ECCO_20th'
modelshort = 'E2'
root = '/home/evankielley/IceDEF/WagnerModel'
condloc = root + '/conditions/' + modelfull + '/'
outloc = root + '/output/' + modelfull + '/'
modelloc = root + '/Model/'
pyOutloc = '/home/evankielley/IceDEF/PyWagnerModel/'

msk = sio.loadmat(condloc + 'mask.mat'); msk = msk['msk'] 
vel = sio.loadmat(condloc + modelshort + '_vels_1992.mat'); vel = vel['vel']
sst = sio.loadmat(condloc + modelshort +'_sst_1992.mat'); sst = sst['sst']
bergdims = sio.loadmat(modelloc + 'bergdims.mat'); bergdims = bergdims['bergdims']; bergdims = np.asarray(bergdims).astype(float)
Laurent_Seed = sio.loadmat(modelloc + 'Laurent_Seed.mat')
Seed_X = Laurent_Seed['Seed_X']; Seed_Y = Laurent_Seed['Seed_Y']
seed_X = np.matlib.repmat(Seed_X, 1, 100); seed_Y = np.matlib.repmat(Seed_Y, 1, 100)
seed_X = seed_X.transpose().flatten(); seed_Y = seed_Y.transpose().flatten()
LAT = vel['latw'] * 1.0; LAT = LAT[0,0]; LAT = np.ravel(LAT); LAT = np.asarray([float(i) for i in LAT])
LON = vel['lonw'] * 1.0; LON = LON[0,0]; LON = np.ravel(LON); LON = np.asarray([float(i) for i in LON])
uwF = vel['uw']; uwF = uwF[0,0]; uwF = np.asarray(uwF).astype(float) 
vwF = vel['vw']; vwF = vwF[0,0]; vwF = np.asarray(vwF).astype(float)
uaF = vel['ua']; uaF = uaF[0,0]; uaF = np.asarray(uaF).astype(float)
vaF = vel['va']; vaF = vaF[0,0]; vaF = np.asarray(vaF).astype(float)
sst = sst[:,:,t]; sst = np.asarray(sst).astype(float)

inFile = outloc + 'output_full.mat'                                                                        
mXIL = sio.loadmat(inFile)['XIL']; mXIL = np.asarray(mXIL).astype(float)                                   
mYIL = sio.loadmat(inFile)['YIL']; mYIL = np.asarray(mYIL).astype(float)                                    
mXI = sio.loadmat(inFile)['mXI']; mXI = np.asarray(mXI).astype(float)                                    
mYI = sio.loadmat(inFile)['mYI']; mYI = np.asarray(mYI).astype(float)                                    
mVOL = sio.loadmat(inFile)['VOL']; mVOL = np.asarray(mVOL).astype(float)                                    
mDVOL = sio.loadmat(inFile)['DVOL']; mDVOL = np.asarray(mDVOL).astype(float)                                  
mL = sio.loadmat(inFile)['mL']; mL = np.asarray(mL).astype(float)                                    
mW = sio.loadmat(inFile)['mW']; mW = np.asarray(mW).astype(float)                                    
mH = sio.loadmat(inFile)['mH']; mH = np.asarray(mH).astype(float)                                    
mUI = sio.loadmat(inFile)['UI']; mUI = np.asarray(mUI).astype(float)                                      
mVI = sio.loadmat(inFile)['VI']; mVI = np.asarray(mVI).astype(float)                                      
mUA = sio.loadmat(inFile)['UA']; mUA = np.asarray(mUA).astype(float)                                      
mVA = sio.loadmat(inFile)['VA']; mVA = np.asarray(mVA).astype(float)                                      
mUW = sio.loadmat(inFile)['UW']; mUW = np.asarray(mUW).astype(float)                                      
mVW = sio.loadmat(inFile)['VW']; mVW = np.asarray(mVW).astype(float)                                      
mTE = sio.loadmat(inFile)['TE']; mTE = np.asarray(mTE).astype(float)                                      
mMemat = sio.loadmat(inFile)['Memat']; mMemat = np.asarray(mMemat).astype(float)                                
mMvmat = sio.loadmat(inFile)['Mvmat']; mMvmat = np.asarray(mMvmat).astype(float)                                
mMbmat = sio.loadmat(inFile)['Mbmat']; mMbmat = np.asarray(mMbmat).astype(float)                                
mLAT = sio.loadmat(inFile)['mLAT']; mLAT = np.asarray(mLAT).astype(float)                                    
mLON = sio.loadmat(inFile)['mLON']; mLON = np.asarray(mLON).astype(float)                                    
mmsk = sio.loadmat(inFile)['mmsk']; mmsk = np.asarray(mmsk).astype(float)                                    
                                                                        

fixedInFile = outloc + 'fixed.mat'
mts = sio.loadmat(fixedInFile)['ts_all']#; mts = np.asarray(mts).astype(float)                                 
mrandoX = sio.loadmat(fixedInFile)['randoX_all']#; mrandoX = np.asarray(mrandoX).astype(float)                        
mrandoY = sio.loadmat(fixedInFile)['randoY_all']#; mrandoY = np.asarray(mrandoY).astype(float)                        

