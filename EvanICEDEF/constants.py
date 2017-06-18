import scipy.io as sio
import numpy as np
import numpy.matlib

R = 6378 * 1e3  ## earth radius in m
rhow = 1027  ## density of water (kg/m^3)
rhoa = 1.2  ## density of air   (kg/m^3)
rhoi = 850  ## density of ice   (kg/m^3)
drho = rhow - rhoi
Cw = 0.9  ## bulk coefficient water  (Bigg et al 1997)
Ca = 1.3  ## bulk coefficient air    (Bigg et al 1997)
om = 7.2921e-5  ## rotation rate of earth (rad/s)
g = np.sqrt(rhoa*drho/rhow/rhoi*(Ca/Cw))  ## gamma = np.sqrt(ca/cw)

Ti0 = -4
Cs1 = 1.5; Cs2 = 0.5; Cs3 = 0.1
CMv1 = 7.62e-3; CMv2 = 1.29e-3; CMe1 = 0.5
CMb1 = 0.58; CMb2 = 0.8; CMb3 = 0.2

bvec = range(8,11)
trajnum = 25            # total number of iceberg trajectories to compute
final_t = 122           # number of input field time steps
startrange = final_t / 2  # input field start range
tres = 3                # time resoln such that "model Dt"="input DT"/tres
DT = 3                  # Input fields time step
Dt = DT / tres            # model timestep in days
dt = Dt * 24 * 3600         # model timestep in seconds
R = 6378 * 1e3            # earth radius in m
dtR = dt/R * 180/np.pi       # need this ratio for distances in "drifting.m"
t = range(0, final_t)      # how long is the run
nt = len(t)*tres            # number of model timesteps
tt = np.linspace(0, len(t)-1,nt)  # model time

modelfull = 'ECCO_20th'
modelshort = 'E2'
root = '/home/evankielley/WagnerModel'
condloc = root + '/conditions/' + modelfull + '/'
outloc = root + '/output/' + modelfull + '/'
modelloc = root + '/Model/'
pyOutloc = '/home/evankielley/IceDEF/EvanICEDEF/'

msk = sio.loadmat(condloc + 'mask.mat'); msk = msk['msk'] 
vel = sio.loadmat(condloc + modelshort + '_vels_1992.mat'); vel = vel['vel']
sst = sio.loadmat(condloc + modelshort +'_sst_1992.mat'); sst = sst['sst']
bergdims = sio.loadmat(modelloc + 'bergdims.mat'); bergdims = bergdims['bergdims']
Laurent_Seed = sio.loadmat(modelloc + 'Laurent_Seed.mat')
Seed_X = Laurent_Seed['Seed_X']; Seed_Y = Laurent_Seed['Seed_Y']
seed_X = np.matlib.repmat(Seed_X, 1, 100); seed_Y = np.matlib.repmat(Seed_Y, 1, 100)
seed_X = seed_X.transpose().flatten(); seed_Y = seed_Y.transpose().flatten()
LAT = vel['latw'] * 1.0; LAT = LAT[0,0]; LAT = np.ravel(LAT)
LON = vel['lonw'] * 1.0; LON = LON[0,0]; LON = np.ravel(LON)
uwF = vel['uw']; vwF = vel['vw']
uwF = uwF[0,0]; vwF = vwF[0,0]
uaF = vel['ua']; vaF = vel['va']  
uaF = uaF[0,0]; vaF = vaF[0,0]
sst = sst[:,:,t]

inFile = outloc + 'output_full.mat'                                                                        
mXIL = sio.loadmat(inFile)['XIL']                                   
mYIL = sio.loadmat(inFile)['YIL']                                   
mVOL = sio.loadmat(inFile)['VOL']                                   
mDVOL = sio.loadmat(inFile)['DVOL']                                 
mUI = sio.loadmat(inFile)['UI']                                     
mVI = sio.loadmat(inFile)['VI']                                     
mUA = sio.loadmat(inFile)['UA']                                     
mVA = sio.loadmat(inFile)['VA']                                     
mUW = sio.loadmat(inFile)['UW']                                     
mVW = sio.loadmat(inFile)['VW']                                     
mTE = sio.loadmat(inFile)['TE']                                     
mMemat = sio.loadmat(inFile)['Memat']                               
mMvmat = sio.loadmat(inFile)['Mvmat']                               
mMbmat = sio.loadmat(inFile)['Mbmat']                               
                                                                        
fixedInFile = outloc + 'fixed.mat'
mts = sio.loadmat(fixedInFile)['ts_all']                                 
mrandoX = sio.loadmat(fixedInFile)['randoX_all']                         
mrandoY = sio.loadmat(fixedInFile)['randoY_all']                         

