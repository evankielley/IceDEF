"""Runs drift simulation using iceberg, metocean, and drift modules from the icedef package.
"""

import pandas as pd
import numpy as np
from datetime import timedelta

class DriftSimulation():
    """This class simulates the drift of an iceberg.

    Notes:
        The ocean and atm objects should contain enough data for the
        intended time and space ranges of the drift simulation.

    Args:
        berg (icedef.iceberg.Iceberg): iceberg object
        ocean (icedef.metocean.ECMWFOcean): ocean object
        atm (icedef.metocean.ECMWFAtm): atmosphere object
        drift (icedef.turnbull.drift): drift function
    """

    OM = 7.2921e-5  # angular velocity of Earth (rad/s)
    RE = 6378e3  # radius of Earth (m)
    RHOA = 1.225 # density of air (kg/m^3)
    RHOW = 1027.5  # density of water (kg/m^3)
    RHOI = 900  # density of iceberg (kg/m^3)

    def __init__(self, berg, ocean, atm, drift):

        self.berg = berg
        self.ocean = ocean
        self.atm = atm
        self.drift = drift

        self.history = None

    @property
    def x_min(self):
        return min(self.ocean.lons[0], self.atm.lons[0])

    @property
    def x_max(self):
        return max(self.ocean.lons[-1], self.atm.lons[-1])

    @property
    def y_min(self):
        return min(self.ocean.lats[0], self.atm.lats[0])

    @property
    def y_max(self):
        return max(self.ocean.lats[-1], self.atm.lats[-1])


    def interpolate(self, t, x, y):

        vcx, vcy = self.ocean.interpolate(t, x, y)
        vax, vay = self.atm.interpolate(t, x, y)

        return vcx, vcy, vax, vay


    def in_bounds(self, x, y):

        if not self.x_min < x < self.x_max:
            print('Iceberg out of bounds')
            return False

        elif not self.y_min < y < self.y_max:
            print('Iceberg out of bounds')
            return False

        else:
            return True

    def meters_to_degrees(self, x, y, dx, dy):
        # x and y are in degrees and dx and dy are in meters

        y1 = y + dy*(180/(np.pi * self.RE))
        x1 = x + dx/(np.cos((((y + y1)/2)*np.pi)/180))*(180/(np.pi*self.RE))

        return x1, y1

    def setup_timestepper(self, nt):

        t = [None]*(nt+1)
        x = np.zeros(nt+1); y = np.zeros(nt+1)
        vx = np.zeros(nt+1); vy = np.zeros(nt+1)
        ax = np.zeros(nt+1); ay = np.zeros(nt+1)

        t[0] = self.berg.t
        x[0] = self.berg.x; y[0] = self.berg.y
        vx[0] = self.berg.vx; vy[0] = self.berg.vy
        ax[0] = self.berg.ax; ay[0] = self.berg.ay

        vcx, vcy = self.ocean.interpolate(t[0], x[0], y[0])
        vax, vay = self.atm.interpolate(t[0], x[0], y[0])

        C = [
            [vcx, vcy, vax, vay],
            [self.berg.M,  self.berg.Ak, self.berg.Ab, self.berg.As, self.berg.At],
            [self.berg.Cdw, self.berg.Cda, self.berg.Csdw, self.berg.Csda],
            [self.OM, self.RHOW, self.RHOA, self.RHOI]
            ]



        return t, x, y, vx, vy, ax, ay, C


    def euler(self, dt, nt):

        t, x, y, vx, vy, ax, ay, C = self.setup_timestepper(nt)

        for i in range(nt):

            ax[i], ay[i] = self.drift(t[i], x[i], y[i], vx[i], vy[i], C)
            vx[i+1] = vx[i] + dt*ax[i]; vy[i+1] = vy[i] + dt*ay[i]
            dx = vx[i+1]*dt; dy = vy[i+1]*dt
            x[i+1], y[i+1] = self.meters_to_degrees(x[i], y[i], dx, dy)
            t[i+1] = t[i] + timedelta(seconds=dt)
            vcx, vcy = self.ocean.interpolate(t[i+1], x[i+1], y[i+1])
            vax, vay = self.atm.interpolate(t[i+1], x[i+1], y[i+1])
            C[0] = [vcx, vcy, vax, vay]

        self.update_history(t, x, y, vx, vy, ax, ay)


    def rk2(self, dt, nt):

        t, x, y, vx, vy, ax, ay, C = self.setup_timestepper(nt)

        for i in range(nt):

            # Stage 1

            ax1, ay1 = self.drift(t[i], x[i], y[i], vx[i], vy[i], C)

            vx1 = vx[i] + 0.5*dt*ax1
            vy1 = vy[i] + 0.5*dt*ay1
            dx1 = dt*vx1
            dy1 = dt*vy1
            x1, y1 = self.meters_to_degrees(x[i], y[i], dx1, dy1)
            t1 = t[i] + 0.5*timedelta(seconds=dt)


            vcx1, vcy1 = self.ocean.interpolate(t1, x1, y1)
            vax1, vay1 = self.atm.interpolate(t1, x1, y1)

            C[0] = [vcx1, vcy1, vax1, vay1]


            # Stage 2

            ax[i], ay[i] = self.drift(t1, x1, y1, vx1, vy1, C)

            vx[i+1] = vx[i] + dt*ax[i]
            vy[i+1] = vy[i] + dt*ay[i]
            dx = dt*vx[i+1]
            dy = dt*vy[i+1]
            x[i+1], y[i+1] = self.meters_to_degrees(x[i], y[i], dx, dy)
            t[i+1] = t[i] + timedelta(seconds=dt)


            vcx, vcy = self.ocean.interpolate(t[i+1], x[i+1], y[i+1])
            vax, vay = self.atm.interpolate(t[i+1], x[i+1], y[i+1])
            C[0] = [vcx, vcy, vax, vay]


        self.update_history(t, x, y, vx, vy, ax, ay)



    def ab2(self, dt, nt):

        t, x, y, vx, vy, ax, ay, C = self.setup_timestepper(nt)

        # Euler

        ax[0], ay[0] = self.drift(t[0], x[0], y[0], vx[0], vy[0], C)
        vx[1] = vx[0] + dt*ax[0]; vy[1] = vy[0] + dt*ay[0]
        dx = vx[1]*dt; dy = vy[1]*dt
        x[1], y[1] = self.meters_to_degrees(x[0], y[0], dx, dy)
        t[1] = t[0] + timedelta(seconds=dt)

        vcx, vcy = self.ocean.interpolate(t[1], x[1], y[1])
        vax, vay = self.atm.interpolate(t[1], x[1], y[1])
        C[0] = [vcx, vcy, vax, vay]


        # Adams-Bashforth 2nd order

        for i in range(1, nt):

            ax[i], ay[i] = self.drift(t[i], x[i], y[i], vx[i], vy[i], C)
            vx[i+1] = vx[i] + dt*(1.5*ax[i] - 0.5*ax[i-1])
            vy[i+1] = vy[i] + dt*(1.5*ay[i] - 0.5*ay[i-1])

            dx = dt*vx[i+1]; dy = dt*vy[i+1]
            x[i+1], y[i+1] = self.meters_to_degrees(x[i], y[i], dx, dy)
            t[i+1] = t[i] + timedelta(seconds=dt)

            vcx, vcy = self.ocean.interpolate(t[i+1], x[i+1], y[i+1])
            vax, vay = self.atm.interpolate(t[i+1], x[i+1], y[i+1])
            C[0] = [vcx, vcy, vax, vay]


        self.update_history(t, x, y, vx, vy, ax, ay)



    def update_history(self, t, x, y, vx, vy, ax, ay):
        self.history = pd.DataFrame({'t': pd.Series(t, dtype='object'),
                                     'x': x, 'y': y, 'vx': vx, 'vy': vy,
                                     'ax': ax, 'ay': ay})


    def trim_list(self, i, *args):
        new_args = []
        for arg in args:
            new_args.append(arg[:i+1])
        return new_args
