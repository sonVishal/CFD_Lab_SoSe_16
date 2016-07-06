#ifndef _INITLB_H_
#define _INITLB_H_
#include "LBDefinitions.h"
#include "helper.h"
#include "computeCellValues.h"
#include <mpi/mpi.h>
#include <stdlib.h>
#include <math.h>

// Read the number of components
void readNumComp(int argc, char *argv[]);


int readParameters(
    int *xlength,                       // Reads domain size. Parameter name: "xlength"
    t_component *c,                     // Contains information on the components of the simulation
    double G[numComp][numComp],         // Matrix giving interaction. name: "G00, G01, G10, G11,..."
    int *procsPerAxis,                  // Number of sub divisions in x,y,z directions
    int *timesteps,                     // Number of timesteps. Parameter name: "timesteps"
    int *timestepsPerPlotting,          // Timesteps between subsequent VTK plots. Parameter name: "vtkoutput"
    int argc,                           // Number of arguments. Should equal 2 (program + name of config file
    char *argv[]                        // argv[1] shall contain the path to the config file
);

// Initialize the fields for all the components
void initialiseProblem(t_component *c, int *flagField, const t_procData * const procData);

// Initialises fields for two components
void initialiseComponents(t_component *c, const t_procData * const procData);

// Initialises fields for one component assuming multiphase.
void initialiseFields(t_component * c, const t_procData * const procData);

#endif
