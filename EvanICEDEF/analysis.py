"""
This module does analyis on the output of an Iceberg simulation.
"""


import matplotlib.pyplot as plt                                         
import pickle                                                           
import scipy.io as sio                                                  
from load_objects import *                                              
                                                                        
def compare_outputs(mXIL,mYIL,bb,pyOutloc,trajnum,nt):                  
    """
    This function compares the output of two separate models.

    Args: 
        mXIL (list): Matlab iceberg x location matrix
        mYIL (list): Matlab iceberg y location matrix
        bb (int): Iceberg size class
        pyOutloc (str): Path to Python model output
        trajnum (int): Number of trajectories per iceberg size class
        nt (int): Number of timesteps per run


    Returns:
        dXIL (list): Difference between X locations between model outputs.
        dYIL (list): Difference between Y locations between model outputs.
    """
    pyXIL, pyYIL = load_objects(pyOutloc,trajnum,nt)                    
    dXIL = mXIL[bb-1,:,:] - pyXIL[:,:]                                  
    dYIL = mYIL[bb-1,:,:] - pyYIL[:,:]                                  
    return dXIL, dYIL 

def plot_model_diff(xDiff,yDiff):
    """
    This function plot the X and Y location difference between models.

    Args:
        xDiff (list): A list of all the X location differences between models.
        yDiff (list): A list od all the Y location differences between models.

    Generates:
        A plot with two subplots; one for x and y differences; respectively.
    """
    plt.subplot(121)
    plt.plot(xDiff.transpose())
    plt.ylabel('Difference in X Location')
    plt.xlabel('Timestep')
    plt.title('Comparing Matlab and Python Models')
    plt.subplot(122)
    plt.plot(yDiff.transpose())
    plt.ylabel('Difference in Y Location')
    plt.xlabel('Timestep')
    plt.title('Comparing Matlab and Python Models')
    plt.show()

def plot_berg_location(xil,yil):
    """
    This function plots the xy position of each iceberg over time.

    Args:
        xil (list): X location of the iceberg.
        yil (list): Y location of the iceberg.

    Generates:
        A plot showing the xy position of icebergs over time.
    """
    plt.plot(xil,yil)                                                       
    plt.xlabel('X Location')                                                   
    plt.ylabel('Y Location')                                                   
    plt.title('Iceberg Drift')                                          
    plt.show() 
