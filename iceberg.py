import random

class Iceberg():                                     
    """                                                                 
    Function for getting mass of 'blocky' (cuboid) icebergs.            
    """                                                                 
    length = 72 + random.gauss(0,5)
    width = 72 + random.gauss(0,5)
    height = 22 + random.gauss(0,5)
    densityIce = 920  # kg/m^3                                          
    mass = densityIce * length * width * height  # kilograms 
