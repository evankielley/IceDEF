import numpy as np
import scipy.io as sio
import numpy.matlib
import numpy.linalg

# Constants
R = 6378*1e3
om = 7.2921e-5
rhow = 1027
rhoa = 1.2
rhoi = 850
drho = rhow - rhoi
Cw = 0.9
Ca = 1.3
gam = np.sqrt(rhoa*drho/rhow/rhoi*(Ca/Cw))

sst0 = -4
Cs1 = 1.5; Cs2 = 0.5; Cs3 = 0.1
CMv1 = 7.62e-3; CMv2 = 1.29e-3; CMe1 = 0.5
CMb1 = 0.58; CMb2 = 0.8; CMb3 = 0.2


# Timesteps
t_final = 122           # number of input field time steps
t_arr = np.arange(t_init, t_final, t_inc)                   # how long is the run
tres = 3                # time resoln such that "model Dt"="input DT"/tres
DT = 3                  # Input fields time step
Dt = DT/tres            # model timestep in days
dt = Dt*24*3600         # model timestep in seconds
nt = len(t)*tres        # number of model timesteps


# Paths
root = '/home/evankielley/IceDEF/NewGeneralModel/'
path2inputs = root + 'Inputs/'
path2vels = path2inputs + 'E2_vels_1992.mat'
path2sst = path2inputs + 'E2_sst_1992.mat'
path2mask = path2inputs + 'mask.mat'


# Input fields
msk = sio.loadmat(path2mask)['msk']
vel = sio.loadmat(path2vels)['vel']
sst = sio.loadmat(path2sst)['sst']
uwF = vel['uw']; uwF = uwF[0,0] 
vwF = vel['vw']; vwF = vwF[0,0]
uaF = vel['ua']; uaF = uaF[0,0]
vaF = vel['va']; vaF = vaF[0,0]
sst = np.asarray(sst).astype(float); sst = sst[:,:,:t_final]

LAT = np.ravel(vel['latw'][0,0]); LAT = np.asarray([float(i) for i in LAT])
LON = np.ravel(vel['lonw'][0,0]); LON = np.asarray([float(i) for i in LON])


# Iceberg Inits
init_coords = []  # lat, lon
init_dims = []  # length, width, height


def main():

    x, y = init_coords[0], init_coords[1]
    l, w, h = init_dims[0], init_dims[1], init_dims[2]

    while t < t_max:

        x_new, y_new, l_new, w_new, h_new = F(t, x, y, l, w, h)

        if x_new > max(x_arr) or x_new < min(x_arr) or y_new > max(y_arr) or y_new < min(y_arr):
            # Iceberg has moved out-of-bounds
            break

        elif l_new <= 0 or w_new <= 0 or h_new <= 0:
            # Iceberg has melted
            break

        else:

            if w_new < 0.85*h_new:
                # Rollover
                w_new, h_new = h_new, w_new

            if w_new > l_new:
                w_new, l_new = l_new, w_new

        t += dt

def F(t,x,y,l,w,h):

    # Data

    t_ = find_nearest(t_arr, t)
    x_ = find_nearest(x_arr, x)
    y_ = find_nearest(y_arr, y)   

    vau = vua_field[x_,y_,t_]
    vav = vav_field[x_,y_,t_]
    vwu = vwu_field[x_,y_,t_]
    vwv = vwv_field[x_,y_,t_]
    sst = sst_field[x_,y_,t_]

    # Drifting

    S = np.pi*((l*w)/(l+w))
    ff = 2*om*np.sin((np.abs(y)*np.pi)/180)
    lam = np.sqrt(2)*Cw*(gam*np.linalg.norm(np.array([vau,vav])))/(ff*S)
    
    if lam < 0.1:
        alpha = lam*(lam**4*(lam**4*(lam**4*(-0.0386699020961393*lam**4 + \
            0.055242717280199) - 0.0883883476483184) + \
            0.176776695296637) - 0.707106781186548)
    else:
        alpha = np.multiply(np.divide(np.sqrt(2),np.power(lam, 3)),(1-np.sqrt(1+np.power(lam,4))))
        
    if lam < 0.6:
        beta = lam**3*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*\
            (0.0153268598203613*lam**4 - 0.0151656272365985) + \
            0.0180267866272764) + 0.0219176256311202) - \
            0.0274446790511418) + 0.0357675015202851) - \
            0.0493731785691779) + 0.0745776683282687) - \
            0.132582521472478) + 0.353553390593274)
    else:
        beta = np.real(np.multiply(np.divide(1.,np.power(lam,3.)),cmath.sqrt(np.multiply((4.+np.power(lam,4.)), \
            cmath.sqrt(1.+np.power(lam,4.)))-3.*np.power(lam,4.)-4.)))

    viu = uw-gam*alpha*va+gam*beta*ua
    viv = vw+gam*alpha*ua+gam*beta*va

    y_new = y + (viv*dt)*(180/(np.pi*R))
    x_new = x + (viu*dt)/(np.cos(np.mean(np.array([y,y_new])*np.pi)/180)*(180/(np.pi*R))

    
    # Melting

    Me = CMe1*(Cs1*np.linalg.norm(vau,vav)**Cs2 + Cs3*np.linalg.norm(vau,vav))
    Mv = CMv1*sst + CMv2*sst**2
    Mb = CMb1*np.power(np.sqrt(np.square(viu-vwu)+np.square(viv-vwv)),CMb2)*(sst - sst0)/l**CMb3

    l_new = l - (Mv + Me)*Dt
    w_new = w - (Mv + Me)*Dt
    h_new = h - Mb*Dt

    return x_new, y_new, l_new, w_new, h_new    


if __name__=="__main__":
    main()
