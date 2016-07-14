import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv("scaleDateStrong.dat", delim_whitespace=True, header=None)
data[6] = data[0] * data[1] * data[2]
data[7] = np.nan
data[8] = np.nan

for i in range(1, data.shape[0]):
    data.loc[i, 7] = float(data.loc[0, 3]) / data.loc[i, 3]
    data.loc[i, 8] = float(data.loc[i, 7]) / data.loc[i, 6]

print data
# columns:
# 0=iproc, 1=jproc, 2=kproc, 3-elapsed_time, 4-#cells, 5-MLUPS, 6-#procs, 7-speedup, 8-efficiency

fontsize = 18
markersize = 18
linewidth = 5
color_marker = "k" # see: http://matplotlib.org/api/colors_api.html

matplotlib.rcParams.update({'font.size': fontsize})

plt.figure()
plt.plot(data[6], data[5], 'o', markersize=markersize, color=color_marker)
plt.xlabel("#processors")
plt.ylabel("MLUPS")
plt.savefig("procsVsMLUPS.png")

plt.figure()
plt.plot(data[6], data[3], 'o', markersize=markersize, color=color_marker)
plt.xlabel("#processors")
plt.ylabel("elapsed time")
plt.savefig("procsVsTime.png")

maxProcs=data.tail(1)[6].values

plt.figure()
plt.plot([0, maxProcs], [0, maxProcs], linewidth=linewidth, color="b") # linear speedup line
plt.plot(data[6], data[7], 'o', markersize=markersize, color=color_marker)
plt.xlabel("#processors")
plt.ylabel("speedup")
plt.savefig("procsVsSpeedup.png")

plt.figure()
plt.plot(data[6], data[8], 'o', markersize=markersize, color=color_marker)
plt.xlabel("#processors")
plt.ylabel("efficiency")
plt.savefig("procsVsEfficiency.png")
