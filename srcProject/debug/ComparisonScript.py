
from paraview import simple
from numpy import genfromtxt
import numpy as np
import os
import datetime

path_debug = os.path.dirname(os.path.abspath(__file__)) # the file should be in /debug/!

reader = simple.OpenDataFile(os.path.join(path_debug, "../../srcLB/debug/REFERENCE_SOLUTION.50.vtk"))
readerRef = simple.OpenDataFile(os.path.join(path_debug, "../pv_files/project_combined.50.pvts"))

fileREF = "REF_SOLUTION.csv"
fileCUR = "CURRENT_SOLUTION.csv"

writer = simple.CreateWriter(os.path.join(path_debug,fileREF), readerRef, Precision=16)
writer.WriteAllTimeSteps = 0
writer.FieldAssociation = "Cells"
writer.UpdatePipeline()

writer = simple.CreateWriter(os.path.join(path_debug,fileCUR), reader, Precision=16)
writer.WriteAllTimeSteps = 0
writer.FieldAssociation = "Cells"
writer.UpdatePipeline()

csv_ref = genfromtxt(os.path.join(path_debug, fileREF), delimiter=',', skip_header=1)
csv_MPI = genfromtxt(os.path.join(path_debug, fileCUR), delimiter=',', skip_header=1)


filename = str(datetime.datetime.now())

#some additional information
nrDifferValues = np.sum(csv_ref != csv_MPI)
differVals = np.abs(csv_ref - csv_MPI)[np.abs(csv_ref - csv_MPI) > 0]
if len(differVals):
	maxDiffer = np.max(differVals)
else:
	maxDiffer = np.nan

normDifferAbs = np.linalg.norm(csv_ref - csv_MPI)
normDifferRel = normDifferAbs / np.linalg.norm(csv_MPI)


print "=============================================================\n"
print "                         RESULT\n"
print "=============================================================\n"

print "Max. value that differs: " + str(maxDiffer)
print "Norm absolute difference: " + str(normDifferAbs)
print "Norm relative difference: " + str(normDifferRel)
print "Nr. values that differ: " + str(nrDifferValues)

if normDifferRel <= 1e-13:
	print "TEST SUCCESSFUL!!!"
else:
	print "TEST FAILED!!!"

# np.set_printoptions(threshold=np.nan)
# print csv_ref - csv_MPI
