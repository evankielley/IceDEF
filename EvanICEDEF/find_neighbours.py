def find_berg_grid(LAT,LON,x,y):
    yi2 = []; xi2 = []

    for indice, lon in reversed(list(enumerate(LON))):
        if lon <= x:
            xi2.append(indice)
            break

    for indice, lon in list(enumerate(LON)):
        if lon > x:
            xi2.append(indice)
            break

    for indice, lat in reversed(list(enumerate(LAT))):
        if lat <= y:
            yi2.append(indice)
            break
    
    for indice, lat in list(enumerate(LAT)):
        if lat > y:
            yi2.append(indice)
            break

    return xi2,yi2            
