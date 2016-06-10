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

#Define increment which the number of processes is increased by.
increment = 2;
numProcs  = 4;
numRuns   = 1;#int(math.floor(numProcs/increment));

#Run simulation the given numRun times.
for i in range(0, numRuns):
    s = "mpirun -np "
    s += str(iproc*jproc*kproc)
    s += " ./lbsim scenarios/cavity.dat"
    os.system(s)
    print s;
