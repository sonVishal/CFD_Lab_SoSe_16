#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "helper.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>
#include <mpi/mpi.h>

/** handles the boundaries in our simulation setup */
void treatBoundary(double * collideField, int const * const flagField, const t_procData * const procData);

#endif
