
#NOTE: this script has to be executed with the executable "pvpython" provided by paraview
#      For me it was in the folder [paraview]/bin  

# in terminal execute:
# [PATH_TO_PVPYTHON] ComparisonScript.py
# /home/daniel/Downloads/ParaView-5.0.1-Qt4-OpenGL2-MPI-Linux-64bit/bin/pvpython ComparisonScript.py


# To run a reference solution:

# 1) go to git branch generate_reference_solution (git checkout generate_reference_solution)
# 2) run the current code with parameter: ./lbsim DEBUG_REF.dat
# 3) There should be a new file in /debug/ folder: DEBUG_REFERENCE_SOLUTION.100.vtk
# 4) go back to master (git checkout master)
# 5) run the MPI solution with "DEBUG_MPI_REF.dat"



from paraview import simple
from numpy import genfromtxt
import numpy as np
import os
import datetime

path_debug = os.path.dirname(os.path.abspath(__file__)) # the file should be in /debug/!

reader = simple.OpenDataFile(os.path.join(path_debug, "DEBUG_REFERENCE_SOLUTION.100.vtk"))
readerRef = simple.OpenDataFile(os.path.join(path_debug, "../pv_files/worksheet4.100.pvts"))

writer = simple.CreateWriter(os.path.join(path_debug,"REF_SOLUTION.csv"), reader)
writer.WriteAllTimeSteps = 0
writer.FieldAssociation = "Cells"
writer.UpdatePipeline()


writer = simple.CreateWriter(os.path.join(path_debug,"MPI_SOLUTION.csv"), reader)
writer.WriteAllTimeSteps = 0
writer.FieldAssociation = "Cells"
writer.UpdatePipeline()

csv_ref = genfromtxt(os.path.join(path_debug, "REF_SOLUTION.csv"), delimiter=',', skip_header=1)
csv_MPI = genfromtxt(os.path.join(path_debug,"MPI_SOLUTION.csv"), delimiter=',', skip_header=1)


filename = str(datetime.datetime.now())
if np.all(csv_ref == csv_MPI):
	print "TEST SUCCESSFUL!!!"
	#open(filename+"__TEST_SUCCESSFULL", 'w+')

else:
	print "TEST FAILED!!!"
	#open(filename+"__TEST_FAILED", 'w+')


