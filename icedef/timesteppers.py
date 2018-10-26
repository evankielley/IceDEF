def euler(f, h, x, **kwargs):
    return x + h * f(v, **kwargs)


def rk2(f, h, x, **kwargs):
    k1 = h * f(x, **kwargs)
    k2 = h * f(x + k1 / 2, **kwargs)
    return x + k2


def rk4(f, h, x, **kwargs):
    k1 = h * f(x, **kwargs)
    k2 = h * f(x + k1 / 2, **kwargs)
    k3 = h * f(x + k2 / 2, **kwargs)
    k4 = h * f(x + k3, **kwargs)
    return x + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6