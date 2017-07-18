import pickle
import numpy as np

def load_objects(inFile,trajnum,nt):
    with open(inFile, 'rb') as input:
        XIL  = np.empty([trajnum,nt]); YIL = np.empty([trajnum,nt])
        i=0
        while True:
            try:
                berg0 = pickle.load(input)
                XIL[i,:] = berg0.coords[:,0]; YIL[i,:] = berg0.coords[:,1]
                i += 1
            except EOFError:
                break
    return XIL,YIL
