import numpy as np
from copy import deepcopy
from icedef import tools


def euler(f, h, *args, **kwargs):
    return h * np.array(f(*args, **kwargs))


def rk2(f, h, *k1_args, **kwargs):

    # k1 stage
    k1 = h * np.array(f(*k1_args, **kwargs))

    # k2 stage
    k2_args = np.array(deepcopy(k1_args))
    k2_args[0] += 0.5 * np.timedelta64(int(h), 's')
    k2_args[2] += 0.5 * tools.dy_to_dlat(k1[1])
    k2_args[1] += 0.5 * tools.dx_to_dlon(k1[0], k2_args[2])
    for i in range(3, len(k1_args)):
        k2_args[i] += 0.5 * k1[i-1]
    k2_args = tuple(k2_args)
    k2 = h * np.array(f(*k2_args, **kwargs))

    return k2


def rk4(f, h, *k1_args, **kwargs):

    # k1 stage
    k1 = h * np.array(f(*k1_args, **kwargs))

    # k2 stage
    k2_args = np.array(deepcopy(k1_args))
    k2_args[0] += 0.5 * np.timedelta64(int(h), 's')
    k2_args[2] += 0.5 * tools.dy_to_dlat(k1[1])
    k2_args[1] += 0.5 * tools.dx_to_dlon(k1[0], k2_args[2])
    for i in range(3, len(k1_args)):
        k2_args[i] += 0.5 * k1[i-1]
    k2_args = tuple(k2_args)
    k2 = h * np.array(f(*k2_args, **kwargs))

    # k3 stage
    k3_args = np.array(deepcopy(k1_args))
    k3_args[0] += 0.5 * np.timedelta64(int(h), 's')
    k3_args[2] += 0.5 * tools.dy_to_dlat(k2[1])
    k3_args[1] += 0.5 * tools.dx_to_dlon(k2[0], k3_args[2])
    for i in range(3, len(k1_args)):
        k3_args[i] += 0.5 * k2[i-1]
    k3_args = tuple(k3_args)
    k3 = h * np.array(f(*k3_args, **kwargs))

    # k4 stage
    k4_args = np.array(deepcopy(k1_args))
    k4_args[0] += np.timedelta64(int(h), 's')
    k4_args[2] += tools.dy_to_dlat(k3[1])
    k4_args[1] += tools.dx_to_dlon(k3[0], k4_args[2])
    for i in range(3, len(k1_args)):
        k4_args[i] += k3[i-1]
    k4_args = tuple(k4_args)
    k4 = h * np.array(f(*k4_args, **kwargs))

    return k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6


def am2(f, h, *args, **kwargs):

    i = len(args[0]) - 1

    if i == 0:
        args0 = tuple([arg[0] for arg in args])
        return h * np.array(f(*args0, **kwargs))
    elif i == 1:
        args0 = tuple([arg[0] for arg in args])
        args1 = tuple([arg[1] for arg in args])
        return h * (3 / 2 * np.array(f(*args1, **kwargs)) -
                    1 / 2 * np.array(f(*args0, **kwargs)))
    else:
        args0 = tuple([arg[-3] for arg in args])
        args1 = tuple([arg[-2] for arg in args])
        args2 = tuple([arg[-1] for arg in args])
        return h * (23 / 12 * np.array(f(*args2, **kwargs)) -
                    4 / 3 * np.array(f(*args1, **kwargs)) +
                    5 / 12 * np.array(f(*args0, **kwargs)))
