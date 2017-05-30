def compare_vals(obsName,obs,expName,exp,step,tol=1e-3):
    diff = obs-exp
    if abs(diff)>tol:
        print('step: {}, diff: {}, {} (obs): {}, {} (exp): {}'.format(step,diff,obsName,obs,expName,exp))

def assert_tol(val1,val2,step,tol=1e-3):
    diff = val1-val2
    if abs(diff)>tol:
        raise AssertionError('step: {}, diff: {}'.format(step,diff))
