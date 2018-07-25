"""Iceberg drift model.

The module can simulate the drifting of an iceberg over one timestep.

The equations of motion used come from the following research paper:

    Turnbull, Ian & Fournier, Nicolas & Stolwijk, Michiel & Fosnaes, Tor & Mcgonigal, David. (2014).
    Operational iceberg drift forecasting in Northwest Greenland.
    Cold Regions Science and Technology.
    110. 10.1016/j.coldregions.2014.10.006.

Todo:
    * Extend drift to include n-layer keel

"""

import numpy as np


def drift(t, x, y, vx, vy, constants):

    """Simulates the drift of an iceberg by calculating the forces that act upon it.

    Notes:

        The arg, constants, is a list of lists that contain the necessary constants for this
        drift function.These constants are:

            constants[0] = vcx, vcy, vax, vay
            constants[1] = m,  Ak, Ab, As, At
            constants[2] = Cdw, Cda, Csdw, Csda
            constants[3] = om, rhow, rhoa, rhoi

        Where:
            vax (float): x-component of air velocity (m/s)
            vay (float): y-component of air velocity (m/s)
            vcx (float): x-component of water velocity (m/s)
            vcy (float): y-component of water velocity (m/s)
            m (float): mass of the iceberg (kg)
            Ak (float): area of the keel face of the iceberg (m^2)
            Ab (float): area of the bottom face of the iceberg (m^2)
            As (float): area of the sail face of the iceberg(m^2)
            At (float): area of the top face of the iceberg (m^2)
            Cdw (float): water drag coefficient for the iceberg
            Cda (float): air drag coefficient for the iceberg
            Csdw (float): water skin drag coefficient for the iceberg
            Csda (float): air skin drag coefficient for the iceberg
            om (float): angular velocity of Earth (rad/s)
            rhow (float): seawater density (kg/m^3)
            rhoa (float): air density (kg/m^3)
            rhoi (float): iceberg density (kg/m^3)

    Args:
        t (datetime.datetime): time of the iceberg
        x (float): longitudinal position of the iceberg (degrees)
        y (float): latitudinal position of the iceberg (degrees)
        vx (float): x-component of the iceberg velocity (m/s)
        vy (float): y-component of the iceberg velocity (m/s)
        constants (list of list of float): various constants needed for drift (see notes)


    Returns:
        ax (float): x-component of iceberg acceleration (m/s^2)
        ay (float): y-component of iceberg acceleration (m/s^2)

    """

    vcx, vcy, vax, vay = constants[0]
    m,  Ak, Ab, As, At  = constants[1]
    Cdw, Cda, Csdw, Csda = constants[2]
    om, rhow, rhoa, rhoi = constants[3]

    # Wind force
    Fax = (0.5*rhoa*Cda*As + rhoa*Csda*At)*abs(vax - vx)*(vax - vx)
    Fay = (0.5*rhoa*Cda*As + rhoa*Csda*At)*abs(vay - vy)*(vay - vy)

    # Water force
    Fwx = (0.5*rhow*Cdw*Ak  + rhow*Csdw*Ab)*abs(vcx - vx)*(vcx - vx)
    Fwy = (0.5*rhow*Cdw*Ak + rhow*Csdw*Ab)*abs(vcy - vy)*(vcy - vy)

    # Coriolis force
    f = 2*om*np.sin(np.deg2rad(y))  # Coriolis parameter
    Fcx = +f*vy*m
    Fcy = -f*vx*m

    # Water pressure gradient force
    vwmx = vcx  # x-component of mean water current down to the iceberg keel (m/s)
    vwmy = vcy  # y-component of mean water current down to the iceberg keel (m/s)
    amwx = 0  # x-component of acceleration (time-derivative) of Vmw (m/s^2)
    amwy = 0  # y-component of acceleration (time-derivative) of Vmw (m/s^2)
    Fwpx = m*(amwx + f*vwmx)
    Fwpy = m*(amwy - f*vwmy)

    # Iceberg acceleration
    ax = (Fax + Fcx + Fwx + Fwpx)/(m + 0.5*m)
    ay = (Fay + Fcy + Fwy + Fwpy)/(m + 0.5*m)

    return ax, ay
