#!/usr/bin/env python

import math
from subprocess import call
import os
from matplotlib import pyplot as plt

###########################################################
#Quick and dirty script to test the scaling of lbsim.     #
# ASSUMES this will not be run with high increment values.#
###########################################################

#TODO:
#   * Add print of speedup
#   * Add print of parallel efficiency
#   * Add second graph with weak speedup
#       * Need to also change xlength in cavity.dat
#   * Save plots to file
#   * Add option to not run C program

#Clear file
f = "scaleData.dat"
g = open(f, "w");
g.close();

#Define starting processes.
iproc    = 1;
jproc    = 1;
kproc    = 1;
xlength0 = 25;
xlength  = 25;
procValue = [iproc, jproc, kproc]

iprocStr = "iProc"
jprocStr = "jProc"
kprocStr = "kProc"
procStr = [iprocStr, jprocStr, kprocStr]

sedStr = ["sed -i '/^" , "/c\\","' scenarios/cavity.dat"]

#Set initial values for procs
for i in range(0,3):
    s = sedStr[0] + procStr[i] + sedStr[1] + procStr[i] + "    " + str(procValue[i]) + sedStr[2]
    os.system(s);

#Set initial values for xlength
s = sedStr[0] + "xlength" + sedStr[1] + "xlength " + "    " + str(xlength) + sedStr[2]
os.system(s);


#Define increment which the number of processes is increased by.

increment = 1;
numRuns   = 2;

#Run simulation the given numRun times.
run = "";
for i in range(0, numRuns):

    #Run program for current proc values
    run = "mpirun -np "
    run += str(procValue[0]*procValue[1]*procValue[2])
    run += " ./lbsim scenarios/cavity.dat"
    os.system(run)

    #Update proc value by increment
    for i in range(0,3):
        procValue[i]+=increment;
        s = sedStr[0] + procStr[i] + sedStr[1] + procStr[i] + "    " + str(procValue[i]) + sedStr[2]
        os.system(s);

    #Update xlength
    xlength = xlength0*procValue[0]*procValue[1]*procValue[2];
    s = sedStr[0] + "xlength" + sedStr[1] + "xlength " + "    " + str(xlength) + sedStr[2]
    os.system(s);


# --------------- #
#    Plotting     #
# --------------- #

# SKIPPING THIS STEP ON CLUSTER

#scaleTest.dat arranged columnwise as follows separated by a space
#iproc jproc kproc elapsedTime numCells MLUPS

#g = open(f)# ;
# line = g.readlines();
# g.close();

# #Holds the variables for each run
# _iproc = []
# _jproc = []
# _kproc = []
# _elapsedTime = []
# _numCells = []
# _MLUPS = []


# #Read each variables from file
# with open(f) as g:
    # for line in g:
        # z = line.split(" ");
        # _iproc.append(int(z[0]));
        # _jproc.append(int(z[1]));
        # _kproc.append(int(z[2]));
        # _elapsedTime.append(float(z[3]));
        # _numCells.append(int(z[4]));
        # _MLUPS.append(float(z[5]));

# numVar = len(_elapsedTime);
# speedup    = [0]*numVar
# efficiency = [0]*numVar



# totProcs = [a*b for a,b in zip(_iproc,_jproc)]
# totProcs = [a*b for a,b in zip(totProcs,_kproc)]

# for i in range(0, len(_elapsedTime)):
    # print(i)
    # speedup[i] = _elapsedTime[0]/ _elapsedTime[i]
    # efficiency[i] = speedup[i]/ totProcs[i]

# plt.plot(totProcs, _MLUPS,"*")
# plt.xlabel('# processes')
# plt.ylabel('MLUPS')
# plt.savefig('procsVSMLUPS.png')

# plt.figure();
# plt.plot(totProcs, _elapsedTime, "*");
# plt.xlabel('# processes')
# plt.ylabel('elapsed time')
# plt.savefig('procsVStime.png')


# plt.figure();
# plt.plot(totProcs, speedup, "*", label='speedup');
# plt.xlabel('# processes')
# plt.ylabel('speedup')
# plt.savefig('procsVSspeedup.png')

# plt.figure();
# plt.plot(totProcs, efficiency, "*", label='efficiency');
# plt.xlabel('# processes')
# plt.ylabel('efficiency')
# plt.savefig('procsVSefficiancy.png')

# plt.show();

