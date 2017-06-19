def compare_vals(obsName,obs,expName,exp,step,tol=1e-3):
    diff = obs-exp
    if abs(diff)>tol:
        print('step: {}, diff: {}, {} (obs): {}, {} (exp): {}'.format(step,diff,obsName,obs,expName,exp))

def assert_tol(val1,val2,step,traj,correct=False,ignore=True,tol=1e-3):
    diff = val1-val2
    if abs(diff)>tol:
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

