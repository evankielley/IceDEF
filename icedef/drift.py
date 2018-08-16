def newtonian_drift(iceberg_velocity, current_velocity, wind_velocity,
                    iceberg_constants):
    """Computes instantaneous iceberg acceleration."""

    # Constants
    Omega = EARTH_ROTATION_RATE
    rhoa = AIR_DENSITY
    rhow = SEAWATER_DENSITY
    Ca = ICEBERG_FORM_DRAG_COEFFICIENT_IN_AIR
    Cw = ICEBERG_FORM_DRAG_COEFFICIENT_IN_WATER
    Cda = ICEBERG_SKIN_DRAG_COEFFICIENT_IN_AIR
    Cdw = ICEBERG_SKIN_DRAG_COEFFICIENT_IN_WATER

    Vwx, Vwy = wind_velocity
    Vcx, Vcy = current_velocity
    # Amwx, Amwy = current_acceleration
    Amwx, Amwy = 0, 0

    As = iceberg_constants['sail_area']
    Ak = iceberg_constants['keel_area']
    At = iceberg_constants['top_area']
    Ab = iceberg_constants['bottom_area']
    M = iceberg_constants['mass']

    phi = iceberg_constants['latitude']

    # Args
    Vx, Vy = iceberg_velocity

    # Wind force
    Fax = (0.5 * rhoa * Ca * As + rhoa * Cda * At) * abs(Vwx - Vx) * (Vwx - Vx)
    Fay = (0.5 * rhoa * Ca * As + rhoa * Cda * At) * abs(Vwy - Vy) * (Vwy - Vy)

    # Current force
    Fwx = (0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * abs(Vcx - Vx) * (Vcx - Vx)
    Fwy = (0.5 * rhow * Cw * Ak + rhow * Cdw * Ab) * abs(Vcy - Vy) * (Vcy - Vy)

    # Coriolis force
    f = 2 * Omega * np.sin(np.deg2rad(phi))
    Fcx = f * M * Vy
    Fcy = -f * M * Vx

    # Water pressure force
    Vmwx = Vcx
    Vmwy = Vcy
    Amwx = Amwx
    Amwy = Amwy
    Fwpx = M * (Amwx + f * Vmwx)
    Fwpy = M * (Amwy - f * Vmwy)

    # Iceberg acceleration
    Ax = (Fax + Fwx + Fcx + Fwpx) / M
    Ay = (Fay + Fwy + Fcy + Fwpy) / M

    return Ax, Ay
