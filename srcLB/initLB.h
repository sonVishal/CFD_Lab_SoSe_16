#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"
#include "LBDefinitions.h"
#include <math.h>


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
    int *xlength,                       /* reads domain size. Parameter name: "xlength" */
    double *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[]                        /* argv[1] shall contain the path to the config file */
);


/* initialises the particle distribution functions and the flagfield */
void initialiseFields(t_component * c, int *flagField, int xlength);

#endif
