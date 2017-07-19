import scipy.io as sio
import numpy as np
import numpy.matlib

# Paths
root = '/home/evankielley/IceDEF/GeneralModel/'
path2inputs = root + 'Inputs/'
path2vels = path2inputs + 'E2_vels_1992.mat'
path2sst = path2inputs + 'E2_sst_1992.mat'
path2bergdims = path2inputs + 'bergdims.mat'
path2mask = path2inputs + 'mask.mat'
path2fixed = path2inputs + 'fixed.mat'
path2seed = path2inputs + 'Laurent_Seed.mat'


# Input fields
msk = sio.loadmat(path2mask)['msk']
vel = sio.loadmat(path2vels)['vel']
sst = sio.loadmat(path2sst)['sst']
#bergdims = sio.loadmat(path2bergdims)['bergdims'] * 1.0
#Laurent_Seed = sio.loadmat(path2seed)
#Seed_X = Laurent_Seed['Seed_X']; Seed_Y = Laurent_Seed['Seed_Y']
#seed_X = np.matlib.repmat(Seed_X, 1, 100); seed_Y = np.matlib.repmat(Seed_Y, 1, 100)
#seed_X = seed_X.transpose().flatten(); seed_Y = seed_Y.transpose().flatten()
LAT = np.ravel(vel['latw'][0,0]); LAT = np.asarray([float(i) for i in LAT])
LON = np.ravel(vel['lonw'][0,0]); LON = np.asarray([float(i) for i in LON])
uwF = vel['uw']; uwF = uwF[0,0] 
vwF = vel['vw']; vwF = vwF[0,0]
uaF = vel['ua']; uaF = uaF[0,0]
vaF = vel['va']; vaF = vaF[0,0]

final_t = 122
sst = np.asarray(sst).astype(float); sst = sst[:,:,:final_t]

# FIXED randomizations (for testing)
mts = sio.loadmat(path2fixed)['ts_all']
mrandoX = sio.loadmat(path2fixed)['randoX_all']
mrandoY = sio.loadmat(path2fixed)['randoY_all']
