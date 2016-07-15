#!/usr/bin/env python

import math
from subprocess import call
import os
#from matplotlib import pyplot as plt

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
iproc = 1;
jproc = 1;
kproc = 1;
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


 #Define increment which the number of processes is increased by.
#   We assume the increment will not be big and we apply a circular 
#   increment for the procs.

increment = 1;
numRuns   = 10;

 #Run simulation the given numRun times.
run = "";
counter = 0;
j = 0; #Index to be changed.
for i in range(0, numRuns):

   # #Circular counter
    if counter == 4:
        counter = 0;
    j = counter%3;

   # #Run program for current proc values
    run = "mpirun -np "
    run += str(procValue[0]*procValue[1]*procValue[2])
    run += " ./lbsim scenarios/cavity.dat"
    os.system(run)

   # #Update proc value
    procValue[j]+=increment
    s = sedStr[0] + procStr[j] + sedStr[1] + procStr[j] + "    " + str(procValue[j]) + sedStr[2]
    os.system(s);
    counter+=1;

