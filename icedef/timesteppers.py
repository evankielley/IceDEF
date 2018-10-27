import numpy as np


def euler(f, h, *args, **kwargs):
    return h * np.array(f(*args, **kwargs))


def rk2(f, h, *args, **kwargs):
    """BROKEN - returns high values"""

    k1 = h * np.array(f(*args, **kwargs))

    args = np.array(args)
    args[0] += np.timedelta64(int(h / 2), 's')

    for i in range(1, len(args)):
        args[i] += k1[i-1] / 2

    args = tuple(args)

    return h * np.array(f(*args, **kwargs))


# def rk4(f, h, x, **kwargs):
#     k1 = h * f(x, **kwargs)
#     k2 = h * f(x + k1 / 2, **kwargs)
#     k3 = h * f(x + k2 / 2, **kwargs)
#     k4 = h * f(x + k3, **kwargs)
#     return x + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6
