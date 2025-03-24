import pandas as pd
from netneurotools import utils, stats
import numpy as np
from scipy.stats import zscore,spearmanr
import os


def corr_spin(x, y, spins, nspins, twotailed=True):

    # convert x and y to arrays to avoid dataframe index bugs
    x = np.array(x)
    y = np.array(y)

    r, p_spearman = spearmanr(x, y)
    null = np.zeros((nspins,))

    if len(x) == spins.shape[0] - 1:  # if insula is missing
        x = np.append(x, np.nan)
        y = np.append(y, np.nan)

    # null correlation
    for i in range(nspins):
        tmp = y[spins[:, i]]
        # convert to dataframe in case insula is missing
        # for better handling of nans
        df = pd.DataFrame(np.stack((x, tmp)).T, columns=['x', 'y'])
        null[i] = df["x"].corr(df["y"])

    p = get_perm_p(r, null, twotailed)

    return r, p_spearman, p, null


def get_perm_p(emp, null, twotailed=True):
    if twotailed:
        return (1 + sum(abs(null - np.nanmean(null)) > abs(emp - np.nanmean(null)))) / (len(null) + 1)
    else:
        return (1 + sum(null > emp)) / (len(null) + 1)


atlas = input('input atlas path:')
coords = utils.get_centroids(atlas, image_space=True)
hemiindex = np.arange(1, 384)
hemiindex = (hemiindex+1) % 2
nspins = 10000

spins033 = stats.gen_spinsamples(coords, hemiid=hemiindex, n_rotate=nspins, seed=1234)


folder_path = input('input folder path:')

for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        x = pd.read_csv(file_path, header=None)
        x = np.array(x)
        x = x[:, 1]
        x = zscore(x)
        for file_name1 in os.listdir(folder_path):
            file_path1 = os.path.join(folder_path, file_name1)
            y = pd.read_csv(file_path1, header=None)
            y = np.array(y)
            y = y[:, 1]
            y = zscore(y)
            info_r, info_p_pearson, info_p, _ = corr_spin(x, y, spins033, nspins, False)
            print(f"{file_name},{file_name1},{info_r},{info_p_pearson},{info_p}\n")
