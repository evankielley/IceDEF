import random

class Iceberg():                                     
    """                                                                 
    Function for getting mass of 'blocky' (cuboid) icebergs.            
    """                                                                 
    sailLength = 72 + random.gauss(0,5)
    sailWidth = 72 + random.gauss(0,5)
    sailHeight = 22 + random.gauss(0,5)
    densityIce = 920  # kg/m^3                                          
    keelWidth = 7
    keelHeight = 110

    mass = densityIce * sailLength * sailWidth * sailHeight  # kilograms 
