import pickle
import scipy.io as sio
from load_objects import *

def compare_outputs(matOutloc,pyOutloc,trajnum,nt):

    matXIL = sio.loadmat(matOutloc + 'E2_B1_full.mat')['XIL']
    matYIL = sio.loadmat(matOutloc + 'E2_B1_full.mat')['YIL']
    print(matXIL.shape)
    pyXIL, pyYIL = load_objects(pyOutloc,trajnum,nt)
    print(pyXIL.shape)
    print(matXIL[:,0])
    print(pyXIL[:,0])
    dXIL = matXIL - pyXIL
    dYIL = matYIL - pyYIL

    return dXIL, dYIL
