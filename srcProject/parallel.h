#ifndef _PARALLEL_H_
#define _PARALLEL_H_
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "helper.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

/* struct that describes the iteration parameters in a certain direction */
typedef struct{
    int startInner, endInner;
    int startOuter, endOuter;
    int fixedValue;
} t_iterPara;

// Wrapper around communicate function to communicate each component
void communicateComponents(double** sendBuffer, double**readBuffer, t_component *c,
                           t_procData const * const procData);

// Performs extract, swap and inject for all the directions
void communicate(double** sendBuffer, double**readBuffer, double* collideField,
                 t_procData const * const procData, int tag);

// Exracts the collide field to the send buffer
void extract( double sendBuffer[], double const * const collideField,
              t_iterPara const * const iterPara, t_procData const * const procData,
              const int direction, int const * const index);

// Swaps the data between neighbors
void swap(double** sendBuffer, double** readBuffer, const t_procData *procData, int direction, int tag);

void swapNoComm(double* collideField, t_iterPara const * const iterPara1,
                t_iterPara const * const iterPara2,	t_procData const * const procData,
                const int direction, int const * const index1, int const * const index2);

// Copys the data from the read buffer to the collide field
void inject(double const * const readBuffer, double* collideField, t_iterPara const * const iterPara,
            t_procData const * const procData, const int direction, int const * const index);

// Density communication functions does the same as for the ones above, but for the density
void communicateDensity(double** sendBuffer, double**readBuffer, double* collideField,
	t_procData const * const procData);

int extractDensity(double sendBuffer[], double const * const rho,
	t_iterPara const * const iterPara, t_procData const * const procData, const int direction);

void swapDensity(double** sendBuffer, double** readBuffer, t_procData const * const procData,
	const int direction, const int bufferSize1, const int bufferSize2);

int injectDensity(double const * const readBuffer, double* rho,
	t_iterPara const * const iterPara, t_procData const * const procData, const int direction);

void swapNoCommDensity(double* rho, t_iterPara const * const iterPara1,
    t_iterPara const * const iterPara2, t_procData const * const procData, const int direction);


// Sets the iteration parameters
void p_setCommIterationParameters(t_iterPara * const iterPara, t_procData const*const procData,
                                  const int direction);

// Assigns the 5 indices that are to be communicated
void p_assignIndices(const int face, int *const index);

// Domain decomposition and setting of neighbors
void domainDecompositionAndNeighbors(t_procData *const procData, const int xlength,
                                     const int * const procsPerAxis);

// Initialise the message passing interface
void initialiseMPI(int *rank, int *numRanks, int argc, char *argv[]);

// Broadcasts the values read by proc 0
void broadcastValues(int rank, int *xlength, t_component *c, double G[numComp][numComp],
                     int *procsPerAxis, int *timesteps, int *timestepsPerPlotting);

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
