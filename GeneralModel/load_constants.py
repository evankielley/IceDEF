import scipy.io as sio
import numpy as np
import numpy.matlib

# Model Constants #####################################################

R = 6378*1e3                                # earth radius in m
rhow = 1027                                 # density of water (kg/m^3)
rhoa = 1.2                                  # density of air (kg/m^3)
rhoi = 850                                  # density of ice (kg/m^3)
drho = rhow - rhoi
Cw = 0.9                    # bulk coefficient water (Bigg et al 1997)
Ca = 1.3                    # bulk coefficient air (Bigg et al 1997)
om = 7.2921e-5              # rotation rate of earth (rad/s)
g = np.sqrt(rhoa*drho/rhow/rhoi*(Ca/Cw))    # gamma = np.sqrt(ca/cw)

Ti0 = -4
Cs1 = 1.5; Cs2 = 0.5; Cs3 = 0.1
CMv1 = 7.62e-3; CMv2 = 1.29e-3; CMe1 = 0.5
CMb1 = 0.58; CMb2 = 0.8; CMb3 = 0.2

bvec = range(1,2) #range(1,11)
trajnum = 1  #25            # total number of iceberg trajectories to compute
final_t = 122           # number of input field time steps
startrange = final_t/2  # input field start range
tres = 3                # time resoln such that "model Dt"="input DT"/tres
DT = 3                  # Input fields time step
Dt = DT/tres            # model timestep in days
dt = Dt*24*3600         # model timestep in seconds
R = 6378*1e3            # earth radius in m
dtR = dt/R*180/np.pi    # need this ratio for distances in drifting

t = range(0, final_t)                   # how long is the run
nt = len(t)*tres                        # number of model timesteps
tt = np.linspace(0, len(t)-1,nt)        # model time
