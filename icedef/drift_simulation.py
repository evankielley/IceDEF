"""This module runs a drift simulation using iceberg and metocean objects and a drift function from the icedef package.
"""

import pandas as pd
import numpy as np
from datetime import timedelta

class DriftSimulation():
    """This class simulates the drift of an iceberg.
    
    Notes:
        The ocean and atm objects should contain enough data for the intended time and space ranges of the drift simulation.
    
    Args:
        berg (icedef.iceberg.Iceberg): iceberg object
        ocean (icedef.metocean.ECMWFOcean): ocean object 
        atm (icedef.metocean.ECMWFAtm): atmosphere object
        drift (icedef.turnbull.drift): drift function
    """
           
    OM = 7.2921e-5  # angular velocity of Earth (rad/s)
    RE = 6378*1e3  # radius of Earth (m)
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

        Vcx, Vcy = self.ocean.interpolate(t, x, y)
        Vax, Vay = self.atm.interpolate(t, x, y)
        
        return Vcx, Vcy, Vax, Vay
    
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
        x = np.empty(nt+1)
        y = np.empty(nt+1)
        vx = np.empty(nt+1)
        vy = np.empty(nt+1)
        vcx = np.empty(nt+1)
        vcy = np.empty(nt+1)
        vax = np.empty(nt+1)
        vay = np.empty(nt+1)
        
        t[0] = self.berg.T
        x[0] = self.berg.X
        y[0] = self.berg.Y
        vx[0] = self.berg.Vx
        vy[0] = self.berg.Vy 
        vcx[0], vcy[0], vax[0], vay[0] = self.interpolate(t[0], x[0], y[0])
        
        constants = [[vcx[0], vcy[0], vax[0], vay[0]],
                     [self.berg.M,  self.berg.Ak, self.berg.Ab, self.berg.As, self.berg.At],
                     [self.berg.Cdw, self.berg.Cda, self.berg.Csdw, self.berg.Csda],
                     [self.OM, self.RHOW, self.RHOA, self.RHOI]] 
        
        return t, x, y, vx, vy, vcx, vcy, vax, vay, constants 
    
    
    def euler(self, dt, nt):
        
        t, x, y, vx, vy, vcx, vcy, vax, vay, constants = self.setup_timestepper(nt)
        
        for i in range(nt):
            
            ax, ay = self.drift(t[i], x[i], y[i], vx[i], vy[i], constants)
            
            vx[i+1] = vx[i] + dt*ax
            vy[i+1] = vy[i] + dt*ay
                        
            dx = vx[i+1]*dt
            dy = vy[i+1]*dt

            # convert position from meters to degrees
            x[i+1], y[i+1] = self.meters_to_degrees(x[i], y[i], dx, dy)
            
            t[i+1] = t[i] + timedelta(seconds=dt)
             
            if not self.in_bounds(x[i+1], y[i+1]):
                t, x, y, vx, vy, vcx, vcy, vax, vay = self.trim_list(i, t, x, y, vx, vy, vcx, vcy, vax, vay)
                break
              
            try: 
                vcx[i+1], vcy[i+1], vax[i+1], vay[i+1] = self.interpolate(t[i+1], x[i+1], y[i+1])
            except ValueError as e:
                print(e)
                t, x, y, vx, vy, vcx, vcy, vax, vay = self.trim_list(i, t, x, y, vx, vy, vcx, vcy, vax, vay)
                break
                  
            constants[0] = [vcx[i+1], vcy[i+1], vax[i+1], vay[i+1]]

        self.update_history(t, x, y, vx, vy, vcx, vcy, vax, vay)
    
            
    def rk2(self, dt, nt):

        t, x, y, vx, vy, vcx, vcy, vax, vay, constants = self.setup_timestepper(nt)
        
        for i in range(nt):
        
            ax, ay = self.drift(t[i], x[i], y[i], vx[i], vy[i], constants)
            half_vx = vx[i] + 0.5*dt*ax 
            half_vy = vy[i] + 0.5*dt*ay
            half_dx = dt*half_vx
            half_dy = dt*half_vy
            half_x, half_y = self.meters_to_degrees(x[i], y[i], half_dx, half_dy)
            half_t = t[i] + 0.5*timedelta(seconds=dt)
         
            if not self.in_bounds(half_x, half_y):
                t, x, y, vx, vy, vcx, vcy, vax, vay = self.trim_list(i, t, x, y, vx, vy, vcx, vcy, vax, vay)
                break
             
            try: 
                half_vcx, half_vcy, half_vax, half_vay = self.interpolate(half_t, half_x, half_y)
            except ValueError as e:
                print(e)
                t, x, y, vx, vy, vcx, vcy, vax, vay = self.trim_list(i, t, x, y, vx, vy, vcx, vcy, vax, vay)
                break


            constants[0] = [half_vcx, half_vcy, half_vax, half_vay]
            half_ax, half_ay = self.drift(half_t, half_x, half_y, half_vx, half_vy, constants)
            
            vx[i+1] = vx[i] + dt*half_ax 
            vy[i+1] = vy[i] + dt*half_ay
            dx = dt*vx[i+1]
            dy = dt*vy[i+1]
            x[i+1], y[i+1] = self.meters_to_degrees(x[i], y[i], dx, dy)
            t[i+1] = t[i] + timedelta(seconds=dt)
         
            if not self.in_bounds(x[i+1], y[i+1]):
                t, x, y, vx, vy, vcx, vcy, vax, vay = self.trim_list(i, t, x, y, vx, vy, vcx, vcy, vax, vay)
                break
            
            try: 
                vcx[i+1], vcy[i+1], vax[i+1], vay[i+1] = self.interpolate(t[i+1], x[i+1], y[i+1])
            except ValueError as e:
                print(e)
                t, x, y, vx, vy, vcx, vcy, vax, vay = self.trim_list(i, t, x, y, vx, vy, vcx, vcy, vax, vay)
                break

            constants[0] = [vcx[i+1], vcy[i+1], vax[i+1], vay[i+1]]
            
        self.update_history(t, x, y, vx, vy, vcx, vcy, vax, vay)
    
    
    def update_history(self, t, x, y, vx, vy, vcx, vcy, vax, vay):  
        self.history = pd.DataFrame({'t': pd.Series(t, dtype='object'), 'x': x, 'y': y, 'vx': vx, 'vy': vy,
                                     'vcx': vcx,'vcy': vcy,'vax': vax,'vay': vay})
        
        
    def trim_list(self, i, *args):
        new_args = []
        for arg in args:
            new_args.append(arg[:i+1])
        return new_args
