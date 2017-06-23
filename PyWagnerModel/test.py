import math
from tabulate import tabulate

def print_var_table(var_list):
    print(tabulate(var_list, headers=['Name', 'Python', 'Matlab', 'Difference']))
    print()  # newline

def compare_vals(obsName,obs,expName,exp,step,tol=1e-3):
    diff = obs-exp
    if abs(diff)>tol:
        print('step: {}, diff: {}, {} (obs): {}, {} (exp): {}'.format(step,diff,obsName,obs,expName,exp))

def assert_tol(name1,val1,name2,val2,step,traj,flag='print',tol=1e-5):
    if math.isnan(val1) or math.isnan(val2):
        print('python: {}, matlab {}, traj: {}, step: {}, diff: {}'.format(val1, val2, traj,step,diff))
        print('This function does not compare NaNs.')
        return val1
    diff = val1-val2
    if abs(diff)>tol:
        if flag=='break':
            raise AssertionError('{}: {}, {}: {}, traj: {}, step: {}, diff: {}'.format(name1,val1,name2,val2,traj,step,diff))
            return val1
        elif flag=='correct':
            val1=val2
            print('{}: {}, {}: {}, traj: {}, step: {}, diff: {}'.format(name1,val1,name2,val2,traj,step,diff))
            return val1
        elif flag=='print':
            print('{}: {}, {}: {}, traj: {}, step: {}, diff: {}'.format(name1,val1,name2,val2,traj,step,diff))
            return val1  
        else:
            print('flag must be one of break, print, or correct')
    else:
        return val1

def assert_tol_matrix(val1,val2,step,traj,correct=False,ignore=False,tol=1e-3):
    diff = val1-val2
    if abs(diff).any()>tol:
        if not correct and not ignore:
            raise AssertionError('traj: {}, step: {}, diff: {}'.format(traj,step,diff))
            return val1
        elif correct and ignore:
            val1=val2
            print('traj: {}, step: {} corrected'.format(traj,step))
            return val1
        else:
            return val1  
    else:
        return val1

