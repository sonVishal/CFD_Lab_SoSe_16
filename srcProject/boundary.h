#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "helper.h"
#include "LBDefinitions.h"
#include "parallel.h"
#include "computeCellValues.h"
#include <stdio.h>
#include <mpi/mpi.h>


typedef struct{
    int x, y, z;
} t_iterParaEdge;

/*Wrapper around tratBoundary to handle every component*/
void treatComponentBoundary(t_component *c,int numComp, int const * const flagField, const t_procData * const procData, double **sendBuffer, double **readBuffer);

/** handles the boundaries in our simulation setup */
void treatBoundary(int const * const flagField, double *collideField, const t_procData * const procData, double **sendBuffer, double **readBuffer);

#endif
