import numpy as np

def find_nearest(array,value):                                          
    idx = (np.abs(array-value)).argmin()                                
    return idx

def find_berg_grid(LAT,LON,x,y):
    #yi2 = []; xi2 = []

    for indice, lon in reversed(list(enumerate(LON))):
        if lon <= x:
            #xi2.append(indice)
            xi2a = indice
            break

    for indice, lon in list(enumerate(LON)):
        if lon > x:
            #xi2.append(indice)
            xi2b = indice
            break

    for indice, lat in reversed(list(enumerate(LAT))):
        if lat <= y:
            #yi2.append(indice)
            yi2a = indice
            break
    
    for indice, lat in list(enumerate(LAT)):
        if lat > y:
            #yi2.append(indice)
            yi2b = indice
            break

    #xi2[1] = xi2[1]+1
    #yi2[1] = yi2[1]+1

    #return xi2,yi2            
    return xi2a,xi2b,yi2a,yi2b
