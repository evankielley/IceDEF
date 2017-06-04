import numpy as np

class Iceberg(object):
    """A chunk of floating ice typically drifting in the ocean.
    
    Args:
        nt (int): the number of timesteps that the iceberg will exist for.

    Attributes:
        dims (float): the dimensions length, width, height, and volume over time. 
        dimsChange (float): the change in dims over time.
        location (float): the location of the iceberg in degrees lon/lat over time.
        melted (bool): a logical that tells whether the iceberg has melted or not.
        outOfBounds (bool): a logical that tells whether the iceberg has drifted out of bounds or not.

    Example:
        berg = Iceberg(nt)
        berg.dims[0,:] = L, W, H, VOL
        berg.dimsChange[0,:] = dL, dW, dH, dVOL
        berg.location[0,:] = lat, lon
        berg.melted = False
        berg.outOfBounds = False  
    """
    def __init__(self, nt):
        """The constructor of the Iceberg class."""
        self.dims = np.empty([nt,4])*np.nan
        self.dimsChange = np.empty([nt,4])*np.nan
        self.location = np.empty([nt,2])*np.nan
        self.melted = False
        self.outOfBounds = False 
