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
        berg.location[0,:] = lon, lat
        berg.melted = False
        berg.outOfBounds = False  
    """
    def __init__(self, nt, dims, coordinates):
        """The constructor of the Iceberg class."""
        nt = int(nt)
        length, width, height, = dims[0], dims[1], dims[2]
        volume = length * width * height
        longitude, latitude = coordinates[0], coordinates[1]
        self.dims = np.multiply(np.empty([nt,4]),np.nan)
        self.dims[0,:] = length, width, height, volume
        self.dimsChange = np.multiply(np.empty([nt,4]),np.nan)
        self.dimsChange[0,:] = 0,0,0,0
        self.location = np.multiply(np.empty([nt,2]),np.nan)
        self.location[0,:] = longitude, latitude
        self.melted = False
        self.outOfBounds = False 
        self.grounded = False
