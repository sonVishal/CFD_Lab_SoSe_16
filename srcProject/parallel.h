#ifndef _PARALLEL_H_
#define _PARALLEL_H_
#include <mpi/mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "helper.h"
#include "LBDefinitions.h"


/* struct that describes the iteration parameters in a certain direction */
typedef struct{
    int startInner, endInner;
    int startOuter, endOuter;
    int fixedValue;
} t_iterPara;

// Wrapper around communicate function to communicate each component
void communicateComponents(double** sendBuffer, double**readBuffer, t_component *c, t_procData const * const procData);

// Performs extract, swap and inject for all the directions
void communicate(double** sendBuffer, double**readBuffer, double* collideField, t_procData const * const procData);

// Exracts the collide field to the send buffer
void extract( double sendBuffer[], double const * const collideField, t_iterPara const * const iterPara, t_procData const * const procData,
              const int direction, int const * const index);

// Swaps the data between neighbors
void swap(double** sendBuffer, double** readBuffer, const t_procData *procData, int direction);

// Copys the data from the read buffer to the collide field
void inject(double const * const readBuffer, double* collideField, t_iterPara *const iterPara, t_procData const * const procData,
            const int direction, int const * const index);

// Sets the iteration parameters
void p_setCommIterationParameters(t_iterPara * const iterPara, t_procData const*const procData, const int direction);

// Assigns the 5 indices that are to be communicated
void p_assignIndices(const int face, int *const index);

// Domain decomposition and setting of neighbors
void domainDecompositionAndNeighbors(t_procData *const procData, const int xlength, const int * const procsPerAxis);

// Initialise the message passing interface
void initialiseMPI(int *rank, int *numRanks, int argc, char *argv[]);

// Broadcasts the values read by proc 0
void broadcastValues(int rank, int *xlength, t_component *c, double G[numComp][numComp],
	double *velocityWall, int *procsPerAxis, int *timesteps, int *timestepsPerPlotting);

// Initialise the buffers
void initialiseBuffers(double *sendBuffer[6], double *readBuffer[6], int const * const xlength,
		int const * const neighbours, int procBufferSize[3]);

// Finalise all the processes and join
void finaliseMPI(t_procData const * const procData);

/* produces a stderr text output  */
void Program_Message(char *txt);

/* produces a stderr textoutput and synchronize all processes */
void Programm_Sync(char *txt);

/* all processes will produce a text output, be synchronized and finished */
void Programm_Stop(char *txt);

#endif
