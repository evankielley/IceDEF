import matplotlib.pyplot as plt                                         
import pickle                                                           
import scipy.io as sio                                                  
from load_objects import *                                              
                                                                        
def compare_outputs(mXIL,mYIL,bb,pyOutloc,trajnum,nt):                  
                                                                        
    pyXIL, pyYIL = load_objects(pyOutloc,trajnum,nt)                    
    print(mXIL.shape)                                                   
    print(pyXIL.shape)                                                  
    dXIL = mXIL[bb-1,:,:] - pyXIL[:,:]                                  
    dYIL = mYIL[bb-1,:,:] - pyYIL[:,:]                                  
                                                                        
    return dXIL, dYIL 

                                                                        
def make_plots(x):                                                      
    plt.plot(x.transpose())                                             
    #plt.plot(x)                                                        
    plt.xlabel('Timestep')                                              
    plt.ylabel('XIL')                                                   
    #plt.ylim([-10,10])                                                 
    plt.title('Iceberg Drift')                                          
    plt.show()                                                          
                                                                        
def make_plots2(x,y):                                                   
    plt.plot(x,y)                                                       
    plt.xlabel('XIL')                                                   
    plt.ylabel('YIL')                                                   
    #plt.ylim([-10,10])                                                 
    plt.title('Iceberg Drift')                                          
    plt.show() 
