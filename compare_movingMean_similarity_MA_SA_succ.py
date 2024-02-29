from operator import index
import numpy as np
import pandas as pd
import os
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
from dtaidistance import dtw
from dtaidistance import dtw_visualisation as dtwvis
import random
 
# Import the data files
succ_MA_csvfile = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results','Moving_mean_succ_MA.csv')
succ_SA_csvfile = os.path.join(r'C:\Users\kamata3\Work\Brain\Brain_network\Stressor_arm\Results','Moving_mean_succ_SA.csv')

succ_MA_df = pd.read_csv(succ_MA_csvfile, header=0)
succ_SA_df = pd.read_csv(succ_SA_csvfile, header=0)

connections = succ_MA_df.columns
for i in range(20):
    v1 = succ_MA_df[connections[i]]
    v2 = succ_SA_df[connections[i]]

    # Similarity Analysis: linear correlation
    rho = np.corrcoef(v1, v2)[0,1]
    #sns.pairplot(data)
    # Similarity Analysis: non-linear correlation
    dist_corr = scipy.spatial.distance.correlation(v1, v2)
    kt = scipy.stats.kendalltau(v1, v2) #  range [-1,1]
    # Similarity Analysis: DTW 
    dtw_dist = dtw.distance(v1, v2)
    path = dtw.warping_path(v1, v2)
    if dtw_dist < 0.4:
        plt.figure(figsize=(8,10))
        dtwvis.plot_warping(v1,v2, path, filename="warp.png")
        random.seed(1)
        for idx in range(len(v2)):
            if random.random() < 0.05:
                v2[idx] += (random.random() - 0.5) / 2
        d, paths = dtw.warping_paths(v1, v2, window=25, psi=2)
        best_path = dtw.best_path(paths)
        dtwvis.plot_warpingpaths(v1, v2, paths, best_path)
        plt.show()
        #plt.close()
        print(f' connection: {connections[i]}, \n Perason_corr: {rho}, distance_corr: {dist_corr} , kt: {kt}, DTW: {dtw_dist}  ')
