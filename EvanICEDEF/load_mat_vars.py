import scipy.io as sio

def load_expected_vars(inFile):
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
    return mXIL,mYIL,mVOL,mDVOL,mUI,mVI,mUA,mVA,mUW,mVW,mTE,mMemat,mMvmat,mMbmat

def load_fixed_vars(inFile):
    mts = sio.loadmat(inFile)['ts_all']
    mrandoX = sio.loadmat(inFile)['randoX_all']
    mrandoY = sio.loadmat(inFile)['randoY_all']
    return mts,mrandoX,mrandoY
