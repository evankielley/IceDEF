import scipy.io as sio
import numpy as np
import numpy.matlib

# Paths ###############################################################
# You MUST tailor these paths to suit the machine that runs that program
root_dir = '/home/evankielley/IceDEF/TestCase/'
input_dir = root_dir + 'Inputs/'
output_dir = root_dir + 'Outputs/'

# Model Constants #####################################################

R = 6378*1e3  ## earth radius in m
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

bvec = range(1,11)
trajnum = 25            # total number of iceberg trajectories to compute
final_t = 122           # number of input field time steps
startrange = final_t/2  # input field start range
tres = 3                # time resoln such that "model Dt"="input DT"/tres
DT = 3                  # Input fields time step
Dt = DT/tres            # model timestep in days
dt = Dt*24*3600         # model timestep in seconds
R = 6378*1e3            # earth radius in m
dtR = dt/R*180/np.pi       # need this ratio for distances in "drifting.m"
t = range(0, final_t)      # how long is the run
nt = len(t)*tres            # number of model timesteps
tt = np.linspace(0, len(t)-1,nt)  # model time


# Read in Input fields ################################################

# Input fields
msk = sio.loadmat(input_dir + 'mask.mat')['msk']
vel = sio.loadmat(input_dir + 'E2_vels_1992.mat')['vel']
sst = sio.loadmat(input_dir +'E2_sst_1992.mat')['sst']
bergdims = sio.loadmat(input_dir + 'bergdims.mat')['bergdims'] * 1.0
Laurent_Seed = sio.loadmat(input_dir + 'Laurent_Seed.mat')
Seed_X = Laurent_Seed['Seed_X']; Seed_Y = Laurent_Seed['Seed_Y']
seed_X = np.matlib.repmat(Seed_X, 1, 100); seed_Y = np.matlib.repmat(Seed_Y, 1, 100)
seed_X = seed_X.transpose().flatten(); seed_Y = seed_Y.transpose().flatten()
LAT = np.ravel(vel['latw'][0,0]); LAT = np.asarray([float(i) for i in LAT])
LON = np.ravel(vel['lonw'][0,0]); LON = np.asarray([float(i) for i in LON])
uwF = vel['uw']; uwF = uwF[0,0] 
vwF = vel['vw']; vwF = vwF[0,0]
uaF = vel['ua']; uaF = uaF[0,0]
vaF = vel['va']; vaF = vaF[0,0]
sst = np.asarray(sst).astype(float); sst = sst[:,:,:final_t] 

# FIXED randomizations (for testing)
fixedInFile = input_dir + 'fixed.mat'
mts = sio.loadmat(fixedInFile)['ts_all']
mrandoX = sio.loadmat(fixedInFile)['randoX_all']
mrandoY = sio.loadmat(fixedInFile)['randoY_all']

# MATLAB outputs (for comparing)
inFile = output_dir + 'WagnerTestCaseOutput.mat'                                                                        
mLAT = sio.loadmat(inFile)['mLAT']
mLON = sio.loadmat(inFile)['mLON']
mmsk = sio.loadmat(inFile)['mmsk']
mXIL = sio.loadmat(inFile)['XIL']
mYIL = sio.loadmat(inFile)['YIL']
mXI = sio.loadmat(inFile)['mXI']
mYI = sio.loadmat(inFile)['mYI']
mVOL = sio.loadmat(inFile)['VOL']
mDVOL = sio.loadmat(inFile)['DVOL']
mL = sio.loadmat(inFile)['mL'] 
mW = sio.loadmat(inFile)['mW']
mH = sio.loadmat(inFile)['mH']
mUa = sio.loadmat(inFile)['mUa']
mUT = sio.loadmat(inFile)['mUT']
mUI = sio.loadmat(inFile)['UI']
mVI = sio.loadmat(inFile)['VI']
mUA = sio.loadmat(inFile)['UA']
mVA = sio.loadmat(inFile)['VA']
mUW = sio.loadmat(inFile)['UW']
mVW = sio.loadmat(inFile)['VW']
mSST = sio.loadmat(inFile)['TE']
mMe = sio.loadmat(inFile)['mMe']
mMv = sio.loadmat(inFile)['mMv']
mMb = sio.loadmat(inFile)['mMb']
mALPHA = sio.loadmat(inFile)['mALPHA']
mBETA = sio.loadmat(inFile)['mBETA']
mGROUNDED = sio.loadmat(inFile)['mGROUNDED']
mOB = sio.loadmat(inFile)['mOB']

