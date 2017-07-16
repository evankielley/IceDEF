import scipy.io as sio
import numpy as np

# Note: 
#   Reading 5 million lines with two integers per line takes about: 
#       46 seconds with numpy.loadtxt, 
#       26 seconds with numpy.genfromtxt, 
#       1 second with pandas.read_csv.
#   Hence we write to csv and later use pandas to read it.

in_filename = 'mask.mat'
in_path = '/home/evankielley/IceDEF/WagnerModel/conditions/ECCO_20th/'

out_filename = 'topography_mask.csv'
out_path = '/home/evankielley/IceDEF/GeneralModel/Inputs/TopographyMask/'

data_type = 'integer'  # must be integer or double

num_fields = 1

if num_fields == 1:
    field = 'msk'
    field = np.asarray(sio.loadmat(in_path + in_filename)[field])

elif num_fields == 2:
    field1 = 'vel'
    field2 = 'uw'
    field = np.asarray(sio.loadmat(in_path + in_filename)[field1][field2])

else:
    print('Error')


if data_type == 'integer':
    np.savetxt(out_path + out_filename, field.astype(int), fmt='%i', delimiter=',')

elif data_type == 'double':
    np.savetxt(out_path + out_filename, field.astype(np.float64), fmt='%f', delimiter=',')

else:
    np.savetxt(out_path + out_filename, field, delimiter=',')

    
