import numpy as np

# Icebergs
num_bergs = 1                              
trajnum = 1                               
init_berg_dims = np.asarray([[600., 500., 400.]]) 
init_berg_coords = np.asarray([[310.,50.]])      


# Time-stepping
final_t = 122           # number of input field time steps
tres = 3                # time resoln such that "model Dt"="input DT"/tres
DT = 3                  # Input fields time step
Dt = DT/tres            # model timestep in days
dt = Dt*24*3600         # model timestep in seconds
R = 6378*1e3            # earth radius in m
dtR = dt/R*180/np.pi    # need this ratio for distances in drifting
t = range(0, final_t)                   # how long is the run
nt = len(t)*tres                        # number of model timesteps
tt = np.linspace(0, len(t)-1,nt)        # model time
