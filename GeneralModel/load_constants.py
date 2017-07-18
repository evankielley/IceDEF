import numpy as np

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

