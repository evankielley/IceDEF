import numpy as np

class Iceberg(object):
    """A chunk of floating ice typically drifting in the ocean.
    
    Args:
        num_steps (inum_steps): the number of timesteps that the iceberg will exist for.

    Attributes:
        dims (float): the dimensions length, width, height, and volume over time. 
        dimsChange (float): the change in dims over time.
        location (float): the location of the iceberg in degrees lon/lat over time.
        melted (bool): a logical that tells whether the iceberg has melted or not.
        outOfBounds (bool): a logical that tells whether the iceberg has drifted out of bounds or not.

    Example:
        berg = Iceberg(num_steps)
        berg.dims[0,:] = L, W, H
        berg.dimsChange[0,:] = dL, dW, dH
        berg.location[0,:] = lon, lat
        berg.melted = False
        berg.outOfBounds = False  
        berg.grounded = False
    """
    def __init__(self, num_steps, init_berg_dims, init_berg_coords):
        """The constructor of the Iceberg class."""
        num_steps = int(num_steps)
        self.dims = np.multiply(np.empty([num_steps,3]),np.nan)
        self.dims[0,:] = init_berg_dims
        self.dimsChange = np.multiply(np.empty([num_steps,3]),np.nan)
        self.dimsChange[0,:] = 0.,0.,0.
        self.coords = np.multiply(np.empty([num_steps,2]),np.nan)
        self.coords[0,:] = init_berg_coords
        self.melted = False
        self.outOfBounds = False 
        self.grounded = False
