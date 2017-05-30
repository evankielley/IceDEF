import numpy as np

class Iceberg(object):
    """
    Info about Iceberg Class
    bergId = 
    dims = 
    dimsChange = 
    location =
    melted =
    outOfBounds = 
    """
    def __init__(self, nt):
        self.dims = np.empty([nt,4])*np.nan
        self.dimsChange = np.empty([nt,4])*np.nan
        self.location = np.empty([nt,2])*np.nan
        self.melted = False
        self.outOfBounds = False 
