import pandas as pd
import urllib
import matplotlib.pyplot as plt

dir_path = 'ftp://data.munroelab.ca/pub/iceberg/beacon/'
dir_contents = urllib.request.urlopen(dir_path).read().splitlines()
filenames = [str(listing.split()[-1])[2:-1] for listing in dir_contents]
csv_filenames = [filename for filename in filenames if filename.startswith('0')
                 and filename.endswith('csv')]
kml_filenames = [filename for filename in filenames if filename.startswith('0')
                 and filename.endswith('kml')]
metadata_filename = filenames[-1]


def get_df(data_dir_path, data_fname):
    df = pd.read_csv(data_dir_path + data_fname)
    df.loc[:, 'DataDate_UTC'] = pd.to_datetime(df['DataDate_UTC'])

    return df


def get_day_idxs(dft, day_inc=5):
    day_idxs = []
    day_js = []
    day_j = 0

    for i in range(len(df)):
        day = (dft[i] - dft[0]).days
        if day == day_j:
            day_idxs.append(i)
            day_js.append(day_j)
            day_j += day_inc

    return day_idxs, day_js


def plot_drift_track(df, data_fname):
    fig = plt.figure()

    x, y = df['Longitude'], df['Latitude']

    plt.scatter(x, y, s=2)

    t = df['DataDate_UTC']
    day_idxs, day_js = get_day_idxs(t)

    for i, day_idx in enumerate(day_idxs):
        plt.text(x[day_idx], y[day_idx], str(day_js[i]), fontsize=16, )

    plt.title(f'Drift Track from {data_fname}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    return fig