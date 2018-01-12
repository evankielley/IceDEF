import numpy as np
import scipy.io as sio
import numpy.matlib
import cmath

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
t0 = 0  # initial time
tn = 122  # total number of timesteps
dt = 1*24*3600  # model timestep in seconds
tf = t0 + tn*dt  # final timestep in seconds  
t_all = np.arange(t0, tn, 1)


# Paths
root = '/home/evankielley/IceDEF/TestCase/'
path2inputs = root + 'Inputs/'
path2vels = path2inputs + 'E2_vels_1992.mat'
path2sst = path2inputs + 'E2_sst_1992.mat'
path2mask = path2inputs + 'mask.mat'
path2outputs = root + 'Outputs/'

# Input fields
msk = sio.loadmat(path2mask)['msk']
vel = sio.loadmat(path2vels)['vel']
sst_field = sio.loadmat(path2sst)['sst']
sst_field = np.asarray(sst_field).astype('Float64'); sst_field = sst_field[:,:,:tn]
vwu_field = np.asarray(vel['uw'][0,0]).astype('Float64')
vwv_field = np.asarray(vel['vw'][0,0]).astype('Float64')
vau_field = np.asarray(vel['ua'][0,0]).astype('Float64')
vav_field = np.asarray(vel['va'][0,0]).astype('Float64')

LAT = np.ravel(vel['latw'][0,0]); LAT = np.asarray([float(i) for i in LAT])
LON = np.ravel(vel['lonw'][0,0]); LON = np.asarray([float(i) for i in LON])


# Iceberg Inits
x_all = LON
y_all = LAT
x0, y0 = 310, 50   # lon, lat
l0, w0, h0 = 600, 500, 400   # length, width, height

def main():

    x, y = x0, y0
    l, w, h = l0, w0, h0
    t = t0
    iceberg = np.array([[t0],[x0],[y0],[l0],[w0],[h0]])

    while t < max(t_all):

        x_new, y_new, l_new, w_new, h_new = iceDEF(t, x, y, l, w, h)

        if x_new > max(x_all) or x_new < min(x_all) or y_new > max(y_all) or y_new < min(y_all):
            # Iceberg out-of-bounds
            print(x_new)
            print(y_new)
            print('out-of-bounds')
            break

        elif l_new <= 0 or w_new <= 0 or h_new <= 0:
            # Iceberg melted
            print('melted')
            break

        else:
            x, y, l, w, h = x_new, y_new, l_new, w_new, h_new
            t += 1
            iceberg_new = np.array([[t],[x],[y],[l],[w],[h]])
            iceberg = np.column_stack((iceberg, iceberg_new))

    save_dict = {
                    't_arr': iceberg[0,:],
                    'x_arr': iceberg[1,:],
                    'y_arr': iceberg[2,:],
                    'l_arr': iceberg[3,:],
                    'w_arr': iceberg[4,:],
                    'h_arr': iceberg[5,:]
                }

    sio.savemat(path2outputs + 'NewPyWagnerTestCaseOutput.mat', save_dict)
    

def find_nearest(array,value):
    """Returns the indice of the closest Lat or Lon to input y or x"""
    value = float(value)
    idx = (abs(array-value)).argmin()
    return idx

def iceDEF(t,x,y,l,w,h):

    print('timestep: {}'.format(t))

    # Data

    t_ = find_nearest(t_all, t)
    x_ = find_nearest(x_all, x)
    y_ = find_nearest(y_all, y)   
    print('t_ = ', t_)
    print('x_ = ', x_)
    print('y_ = ', y_)

    vau = vau_field[x_,y_,t_]
    vav = vav_field[x_,y_,t_]
    vwu = vwu_field[x_,y_,t_]
    vwv = vwv_field[x_,y_,t_]
    sst = sst_field[x_,y_,t_]
    print('vau = ', vau)
    print('vav = ', vav)
    print('vwu = ', vwu)
    print('vwv = ', vwv)
    print('sst = ', sst)

    # Drifting

    S = np.pi*((l*w)/(l+w))
    ff = 2*om*np.sin((np.abs(y)*np.pi)/180)
    lam = np.sqrt(2)*Cw*(gam*np.sqrt(vau**2 + vav**2))/(ff*S)
    print('S = {0:.15f}'.format(S))
    print('ff = {0:.15f}'.format(ff))
    print('lam = {0:.15f}'.format(lam))
    
    if lam < 0.1:
        print('Taylor approx used for alpha')
        alpha = lam*(lam**4*(lam**4*(lam**4*(-0.0386699020961393*lam**4 + \
            0.055242717280199) - 0.0883883476483184) + \
            0.176776695296637) - 0.707106781186548)
    else:
        alpha = np.multiply(np.divide(np.sqrt(2),np.power(lam, 3)),(1-np.sqrt(1+np.power(lam,4))))
        
    if lam < 0.6:
        print('Taylor approx used for beta')
        beta = lam**3*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*(lam**4*\
            (0.0153268598203613*lam**4 - 0.0151656272365985) + \
            0.0180267866272764) + 0.0219176256311202) - \
            0.0274446790511418) + 0.0357675015202851) - \
            0.0493731785691779) + 0.0745776683282687) - \
            0.132582521472478) + 0.353553390593274)
    else:
        beta = np.real(np.multiply(np.divide(1.,np.power(lam,3.)),cmath.sqrt(np.multiply((4.+np.power(lam,4.)), \
            cmath.sqrt(1.+np.power(lam,4.)))-3.*np.power(lam,4.)-4.)))

    print('alpha = {0:.15f}'.format(alpha))
    print('beta = {0:.15f}'.format(beta))

    viu = vwu + gam*(-alpha*vav + beta*vau)
    viv = vwv + gam*(alpha*vau + beta*vav)

    print('viu = {0:.15f}'.format(viu))
    print('viv = {0:.15f}'.format(viv))

    y_new = y + (viv*dt)*(180/(np.pi*R))
    x_new = x + (viu*dt)/(np.cos((((y + y_new)/2)*np.pi)/180))*(180/(np.pi*R))

    print('x_new = {0:.15f}'.format(x_new))
    print('y_new = {0:.15f}'.format(y_new))
    
    # Melting

    Me = CMe1*(Cs1*np.sqrt(vau**2 + vav**2)**Cs2 + Cs3*np.sqrt(vau**2 + vav**2))
    Mv = CMv1*sst + CMv2*sst**2
    Mb = CMb1*np.power(np.sqrt(np.square(viu-vwu)+np.square(viv-vwv)),CMb2)*(sst - sst0)/l**CMb3

    print('Me = {0:.15f}'.format(Me))
    print('Mv = {0:.15f}'.format(Mv))
    print('Mb = {0:.15f}'.format(Mb))

    l_new = l - (Mv + Me)*(dt/(24*3600))  # convert dt from secs to days
    w_new = w - (Mv + Me)*(dt/(24*3600))
    h_new = h - Mb*(dt/(24*3600))

    if w_new < 0.85*h_new:
        # Rollover
        print('rollover')
        w_new, h_new = h_new, w_new

    if w_new > l_new:
        # Ensure l is greater than w
        print('swap l and w')
        w_new, l_new = l_new, w_new

    print('l_new = {0:.15f}'.format(l_new))
    print('w_new = {0:.15f}'.format(w_new))
    print('h_new = {0:.15f}'.format(h_new))

    return x_new, y_new, l_new, w_new, h_new    


if __name__=="__main__":
    main()
