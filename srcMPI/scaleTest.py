import math
from subprocess import call
import os

#Script to test the scaling of lbsim.
# ASSUMES this will not be run with high increment values.

#Clear file
f = open("scaleData.dat", "w");
f.close();

#Define starting processes.
iproc = 2;
jproc = 2;
kproc = 2;
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
#   We assume the increment will not be big and apply a circular 
#   increment for the procs.
increment = 2;
numProcs  = 4;
numRuns   = 1;#int(math.floor(numProcs/increment));

#TODO: Make function that sets procs.
#os.system("sed -i '/^iProc/c\This line is removed by the admin.' scenarios/cavity.dat")
#Run simulation the given numRun times.
counter = 0;
for i in range(0, numRuns):
    s = "mpirun -np "
    s += str(iproc*jproc*kproc)
    s += " ./lbsim scenarios/cavity.dat"
    #os.system(s)
    print s;


def setProcs(iproc, jproc, kproc):
    print iproc
    print jproc
    print kproc
    #TODO: Sed for all procs and set them.
