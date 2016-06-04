#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
    int *xlength,                       /* reads domain size. Parameter name: "xlength" */
    double *tau,                        /* relaxation parameter tau. Parameter name: "tau" */
    double *velocityWall,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
    int *iProc,                         /* Number of sub divisions in x direction*/
    int *jProc,                         /* Number of sub divisions in y direction*/
    int *kProc,                         /* Number of sub divisions in z direction*/
    int *timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
    int *timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
    int argc,                           /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[]                        /* argv[1] shall contain the path to the config file */
);


/* initialises the particle distribution functions and the flagfield */
void initialiseFields(double *collideField, double *streamField,int *flagField, int xlength, int rank, int number_of_ranks);

// Initialise the message passing interface
void initialiseMPI(int *rank, int *number_of_ranks, int argc, char *argv[]);

// Initialise the buffers
void initialiseBuffers(double *sendBuffer[6], double *readBuffer[6], int xlength);

// Finalise all the processes and join
void finaliseMPI();

#endif
