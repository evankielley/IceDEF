import numpy as np

def find_nearest(array,value):                                          
    idx = (np.abs(array-value)).argmin()                                
    return idx

def find_berg_grid(LAT,LON,x,y):

    for indice, lon in reversed(list(enumerate(LON))):
        if lon <= x:
            xi2a = indice
            break

    for indice, lon in list(enumerate(LON)):
        if lon > x:
            xi2b = indice
            break

    for indice, lat in reversed(list(enumerate(LAT))):
        if lat <= y:
            yi2a = indice
            break
    
    for indice, lat in list(enumerate(LAT)):
        if lat > y:
            yi2b = indice
            break

    return xi2a,xi2b,yi2a,yi2b
