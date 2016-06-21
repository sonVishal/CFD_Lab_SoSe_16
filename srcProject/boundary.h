#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "helper.h"
#include "LBDefinitions.h"
#include "parallel.h"
#include "computeCellValues.h"
#include <stdio.h>
#include <mpi/mpi.h>

/** handles the boundaries in our simulation setup */
void treatBoundary(t_component *c, int const * const flagField, double *collideField, const t_procData * const procData, double **sendBuffer, double **readBuffer);

#endif
