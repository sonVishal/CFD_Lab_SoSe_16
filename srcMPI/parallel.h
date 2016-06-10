#ifndef _PARALLEL_H_
#define _PARALLEL_H_
#include <mpi/mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "helper.h"
#include "boundary.h"
#include "LBDefinitions.h"

typedef struct{
    int startInner, endInner;
    int startOuter, endOuter;
    int fixedValue;

} t_iterPara;

// Performs extract, swap and inject for all the directions
void communicate(double** sendBuffer, double**readBuffer, double* collideField, const t_procData *procData);

// Exracts the collide field to the send buffer
void extract( double sendBuffer[], double* collideField, const t_iterPara *iterPara, const t_procData *procData,
              int direction, int* index);

// Swaps the data between neighbors
void swap(double** sendBuffer, double** readBuffer, const t_procData *procData, int direction);

// Copys the data from the read buffer to the collide field
void inject(double readBuffer[], double* collideField, t_iterPara *iterPara, const t_procData *procData,
            int direction, int *index);

// Sets the iteration parameters
void p_setCommIterationParameters(t_iterPara *iterPara, const t_procData *procData, const int direction);

// Assigns the 5 indices that are to be communicated
void p_assignIndices(int direction, int *index);

// Domain decomposition and setting of neighbours
void p_domainDecompositionAndNeighbors(t_procData *procData, const int xlength, const int * const procsPerAxis);

// Initialise the message passing interface
void initialiseMPI(int *rank, int *numRanks, int argc, char *argv[]);

// Initialise the buffers
void initialiseBuffers(double *sendBuffer[6], double *readBuffer[6], int *xlength, int *neighbours,
		int procBufferSize[3]);

// Finalise all the processes and join
void finaliseMPI();

/* produces a stderr text output  */
void Program_Message(char *txt);

/* produces a stderr textoutput and synchronize all processes */
void Programm_Sync(char *txt);

/* all processes will produce a text output, be synchronized and finished */
void Programm_Stop(char *txt);

#endif
