import matplotlib.pyplot as plt
import numpy as np

def plot1(iip_berg, mod_berg):
    tol = 0.1  # tstep
    for i,t in enumerate(iip_berg.t2000):
        diff = abs(np.asarray(mod_berg.t2000) - t)
        for j,diff_t in enumerate(diff):
            if diff_t < tol:
                break
        plt.plot(iip_berg.lats[i], iip_berg.lons[i], marker='o', color='red')
        plt.text(iip_berg.lats[i], iip_berg.lons[i], '{0:.1f}'.format(t-iip_berg.t2000[0]))
        plt.plot(mod_berg.lats[j], mod_berg.lons[j], marker='o', color='yellow')
        plt.text(mod_berg.lats[j], mod_berg.lons[j], '{0:.1f}'.format(t-iip_berg.t2000[0]))
    plt.plot(iip_berg.lats, iip_berg.lons, label='observed', color='red')
    plt.plot(mod_berg.lats, mod_berg.lons, label='computed', color='orange')
    plt.legend()
    plt.xlabel('Latitude'); plt.ylabel('Longitude')
    plt.title('Iceberg: {}\nStart time: {}'.format(iip_berg.id_num, iip_berg.times[0]))
    i = 0
    filename = './drift_track_{}'.format(iip_berg.id_num)
    i = 0
    while True:
        i += 1
        newname = '{}_{:d}.png'.format(filename, i)
        if os.path.exists(newname):
            continue
        plt.savefig(newname)
        break
    plt.show()


def plot_return(iip_berg, mod_berg):
    f = plt.figure()
    tol = 0.1  # tstep
    for i,t in enumerate(iip_berg.t2000):
        diff = abs(np.asarray(mod_berg.t2000) - t)
        for j,diff_t in enumerate(diff):
            if diff_t < tol:
                break
        plt.plot(iip_berg.lons[i], iip_berg.lats[i], marker='o', color='red')
        plt.text(iip_berg.lons[i], iip_berg.lats[i], '{0:.1f}'.format(t-iip_berg.t2000[0]))
        plt.plot(mod_berg.lons[j], mod_berg.lats[j], marker='o', color='yellow')
        plt.text(mod_berg.lons[j], mod_berg.lats[j], '{0:.1f}'.format(t-iip_berg.t2000[0]))
    plt.plot(iip_berg.lons, iip_berg.lats, label='observed', color='red')
    plt.plot(mod_berg.lons, mod_berg.lats, label='computed', color='orange')
    plt.legend()
    plt.xlabel('Longitude'); plt.ylabel('Latitude')
    plt.title('Iceberg: {}\nStart time: {}'.format(iip_berg.id_num, iip_berg.times[0]))
    return f


def plot_return_size_vary(iip_berg, mod_berg_gr, mod_berg_bb, mod_berg_sm, mod_berg_med, 
                         mod_berg_lg, mod_berg_vlg, ind=None):
    f = plt.figure()
    tol = 0.1  # tstep
    for i,t in enumerate(iip_berg.t2000):
        diff = abs(np.asarray(mod_berg_med.t2000) - t)
        for j,diff_t in enumerate(diff):
            if diff_t < tol:
                break
        plt.plot(iip_berg.lons[i], iip_berg.lats[i], marker='o', color='red')
        plt.text(iip_berg.lons[i], iip_berg.lats[i], '{0:.1f}'.format(t-iip_berg.t2000[0]))
        plt.plot(mod_berg_med.lons[j], mod_berg_med.lats[j], marker='o', color='yellow')
        plt.text(mod_berg_med.lons[j], mod_berg_med.lats[j], '{0:.1f}'.format(t-iip_berg.t2000[0]))
    plt.plot(iip_berg.lons, iip_berg.lats, label='observed', color='red')
    plt.plot(mod_berg_gr.lons, mod_berg_gr.lats, label='gr', color='orange')
    plt.plot(mod_berg_bb.lons, mod_berg_bb.lats, label='bb', color='green')
    plt.plot(mod_berg_sm.lons, mod_berg_sm.lats, label='sm', color='blue')
    plt.plot(mod_berg_med.lons, mod_berg_med.lats, label='med', color='black')
    plt.plot(mod_berg_lg.lons, mod_berg_lg.lats, label='lg', color='purple')
    plt.plot(mod_berg_vlg.lons, mod_berg_vlg.lats, label='vlg', color='yellow')
    plt.legend()
    plt.xlabel('Longitude'); plt.ylabel('Latitude')
    if ind is not None:
        plt.title('Index: {}, Iceberg: {}\nStart time: {}'.format(ind,iip_berg.id_num, iip_berg.times[0]))
    else:
        plt.title('Iceberg: {}\nStart time: {}'.format(iip_berg.id_num, iip_berg.times[0]))
    return f




