"""This module tests the likeness between two seperate porgrams"""
import math
from tabulate import tabulate

def print_var_table(var_list,traj,timestep):
    """This function prints a neat table of variables specified in an input array """
    print("trajectory: {}, timestep: {}".format(traj,timestep))
    print(tabulate(var_list, headers=['Name', 'Python', 'Matlab', 'Difference']))
    print()  # newline

def assert_tol(name1,val1,name2,val2,step,traj,flag='print',tol=1e-5):
    """This function uses several different flags (print, break, correct, and ignore) to compare input values"""

    if math.isnan(val1) or math.isnan(val2):
        print(output_str)
        print('This function does not compare NaNs.')
        return val1
    diff = val1-val2
    output_str = 'python: {}, matlab {}, traj: {}, step: {}, diff: {}'.format(val1, val2, traj,step,diff)

    if abs(diff)>tol:
        if flag=='break':
            raise AssertionError(output_str)
            return val1
        elif flag=='correct':
            val1=val2
            print(output_str)
            return val1
        elif flag=='print':
            print(output_str)
            return val1  
        elif flag=='ignore':
            return val1
        else:
            return val1
    else:
        return val1

def assert_tol_matrix(val1,val2,step,traj,flag='ignore',tol=1e-3):
    """This function uses several different flags (print, break, correct, and ignore) to compare input matrices"""
    diff = val1-val2
    output_str = 'traj: {}, step: {}, diff: {}'.format(traj,step,diff)
    if abs(diff).any()>tol:
        if flag=='break':
            raise AssertionError(output_str)
            return val1
        elif flag=='correct':
            val1=val2
            print(output_str)
            return val1
        elif flag=='print':
            print(output_str)
            return val1  
        elif flag=='ignore':
            return val1
        else:
            return val1
    else:
        return val1
