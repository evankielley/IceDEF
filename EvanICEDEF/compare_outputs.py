import pickle
import scipy.io as sio
from load_objects import *

def compare_outputs(mXIL,mYIL,bb,pyOutloc,trajnum,nt):

    pyXIL, pyYIL = load_objects(pyOutloc,trajnum,nt)
    print(mXIL.shape)
    print(pyXIL.shape)
    dXIL = mXIL[bb-1,:,:] - pyXIL[:,:]
    dYIL = mYIL[bb-1,:,:] - pyYIL[:,:]

    return dXIL, dYIL
