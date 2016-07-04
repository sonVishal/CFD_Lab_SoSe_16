#ifndef _INITLB_H_
#define _INITLB_H_
#include "LBDefinitions.h"
#include "helper.h"
#include "computeCellValues.h"
#include <mpi/mpi.h>
#include <stdlib.h>
#include <math.h>

//// Read the number of components
void readNumComp(int argc, char *argv[]);

//[> reads the parameters for the lid driven cavity scenario from a config file <]

int readParameters(
    int *xlength,                       // reads domain size. Parameter name: "xlength"
    t_component *c,
    double G[numComp][numComp],         // Matrix giving interaction. name: "G00, G01, G10, G11,..."
    int *procsPerAxis,                  // Number of sub divisions in x,y,z directions
    int *timesteps,                     // number of timesteps. Parameter name: "timesteps"
    int *timestepsPerPlotting,          // timesteps between subsequent VTK plots. Parameter name: "vtkoutput"
    int argc,                           // number of arguments. Should equal 2 (program + name of config file
    char *argv[]                        // argv[1] shall contain the path to the config file
);

//[> Initialize the fields for all the components <]
void initialiseComponents(t_component *c, int *flagField, const t_procData * const procData);

/* initialises the particle distribution functions and the flagfield */
void initialiseFields(t_component * c, const t_procData * const procData);

#endif
