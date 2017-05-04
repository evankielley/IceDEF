import random
import numpy as np
import matplotlib.pyplot as plt
from iceberg import Iceberg

iceberg = Iceberg()

def coriolis_parameter(ϕ=50):
    """
    calculate the Coriolis parameter f
        f = 2 Ω sin(ϕ)
    for a given latitude ϕ in degrees North
    """
    Ω = 7.2921e-5 # rotation rate of earth in rad/s
    f = 2*Ω*np.sin(np.deg2rad(ϕ)) # calculate the Coriolis parameter
    
    return f

def Fa(iceberg, Vx, Vy, t):
    """
    Calculate the air drag
    """
    
    # air velocity, "wind"
    # nominal value (diurnal? wind for 24 hr cycle)
    #Vax = 0*60.0*np.cos(1*np.pi*t/86400)
    Vax = 5*np.cos((2*np.pi*t)/86400) + 10 + random.gauss(0, 1)  # mean, std
    Vay = 0
    
    # calculate air density
    # TODO: effect of temperature and humidity
    ρa = 1.225 # kg/m^3
    
    Ca = 0.1
    As = iceberg.sailWidth * iceberg.sailHeight # (maximum width (m) * sail height (m))
    
    # TODO: add skin drag
    Fax = (0.5 * ρa * Ca * As) * abs(Vax - Vx) * (Vax - Vx)
    Fay = (0.5 * ρa * Ca * As) * abs(Vay - Vy) * (Vay - Vy)
    
    return Fax, Fay

def Fw(iceberg, Vx, Vy, t):
    """
    Calculate the water drag
    """
    Vcx = 0
    Vcy = 0
    
    ρw = 1027.5 # kg/m^3
    Cw = 2.5 #
    Ak = iceberg.keelWidth * iceberg.keelHeight # width (m) * keel depth (m)
    
    Cdw = 5.0e-4 #
    
    # TODO: extend to n-layer keel model
    # TODO: add skin drag
    Fwx = (0.5 * ρw * Cw * Ak *abs(Vcx - Vx)*(Vcx-Vx)) +(0)
    Fwy = (0.5 * ρw * Cw * Ak *abs(Vcy - Vy)*(Vcy-Vy)) +(0)
    
    return Fwx, Fwy

def Fc(iceberg, Vx, Vy, t):
    """
    Calculate the Coriolis force
    """
    f = coriolis_parameter()
    
    Fcx = + f * Vy * iceberg.mass
    Fcy = - f * Vx * iceberg.mass
    
    return Fcx, Fcy

def Fwp(iceberg, Vx, Vy, t):
    """
    Calculate the water pressure gradient
    """
    # Mean water current down to the iceberg keel
    Vwmx = 0
    Vwmy = 0
    
    # acceleration (time-derivative) of Vmw
    Amwx = 0
    Amwy = 0
    
    f = coriolis_parameter()
    
    Fwpx = iceberg.mass*(Amwx + f*Vwmx)
    Fwpy = iceberg.mass*(Amwy - f*Vwmy)
    
    return Fwpx, Fwpy

def calc_acceleration(iceberg, Vx, Vy, t):
    """
    """
    
    Fax, Fay = Fa(iceberg, Vx, Vy, t)
    Fcx, Fcy = Fc(iceberg, Vx, Vy, t)
    Fwx, Fwy = Fw(iceberg, Vx, Vy, t)
    Fwpx, Fwpy = Fwp(iceberg, Vx, Vy, t)
    
    # added mass
    Ma = 0.5 * iceberg.mass
                
    ax = (Fax + Fcx + Fwx + Fwpx) / (iceberg.mass + Ma)
    ay = (Fay + Fcy + Fwy + Fwpy) / (iceberg.mass + Ma)
    
    return ax, ay

def solve(iceberg, dt=60.0):
    """
    Solve the equations of motion with
        dt time step in seconds
        
        iceberg is an Iceiceberg
    """
    tmax_hours = 72
    tmax = tmax_hours*60*60 # hours to seconds
    
    # total number of timesteps
    N = round(tmax/dt)

    # allocate memory for arrays
    x = np.zeros(N)
    y = np.zeros(N)
    vx = np.zeros(N)
    vy = np.zeros(N)
    ax = np.zeros(N)
    ay = np.zeros(N)
    t = np.zeros(N)

    # initial values
    x[0] = 0
    y[0] = 0
    #vx[0] = 0
    #vy[0] = 0
    vx[0] = 0.1 + random.gauss(0, 0.01)  # mean, std 
    vy[0] = 0.1 + random.gauss(0, 0.01)
    ax[0], ay[0] = calc_acceleration(iceberg, vx[0], vy[0], t[0])
    t[0] = 0

    # integrate numerically
    for i in range(N-1):
        
        t[i+1] = t[i] + dt
            
        if i < 1:
            # explicit Euler forward scheme
            vx[i+1] = vx[i] + dt*ax[i]
            vy[i+1] = vy[i] + dt*ay[i]
        
            ax[i+1], ay[i+1] = calc_acceleration(iceberg, vx[i+1], vy[i+1], t[i+1])
        elif i < 3:
            # second order Adams Bashforth
            vx[i+1] = vx[i] + dt*(1.5*ax[i]-0.5*ax[i-1])
            vy[i+1] = vy[i] + dt*(1.5*ay[i]-0.5*ay[i-1])
            
            ax[i+1], ay[i+1] = calc_acceleration(iceberg, vx[i+1], vy[i+1], t[i+1])
        else:
            # fourth order Adams Bashforth, predictor-corrector
            vx[i+1] = vx[i] + dt/24*(55*ax[i]-59*ax[i-1]+37*ax[i-2]-9*ax[i-3])
            vy[i+1] = vy[i] + dt/24*(55*ay[i]-59*ay[i-1]+37*ay[i-2]-9*ay[i-3])    
            
            ax[i+1], ay[i+1] = calc_acceleration(iceberg, vx[i+1], vy[i+1], t[i+1])
            
            vx[i+1] = vx[i] + dt/24*(9*ax[i+1]+19*ax[i]-5*ax[i-1]+ax[i-2])
            vy[i+1] = vy[i] + dt/24*(9*ay[i+1]+19*ay[i]-5*ay[i-1]+ay[i-2])
            
            ax[i+1], ay[i+1] = calc_acceleration(iceberg, vx[i+1], vy[i+1], t[i+1])
            
        x[i+1] = x[i] + dt*vx[i+1]
        y[i+1] = y[i] + dt*vy[i+1]
  
    return x, y, t

def plot(x, y, label=''):
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.plot(x, y, label=label)

fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis

for i in range(0, 10):
    iceberg = Iceberg()
    x, y, t = solve(iceberg)
    ax.plot(x,y)

plt.title("Varying Iceberg Size, Initial Drift, and Air Speed")
fig.savefig('./driftTrack.png')   # save the figure to file
plt.close(fig)

