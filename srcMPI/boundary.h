#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "helper.h"
#include "debug.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>
#include <mpi/mpi.h>

/** handles the boundaries in our simulation setup */
void treatBoundary(double * collideField, int const * const flagField, const t_procData * const procData);

void p_setIterationParameters(int *endOuter, int *endInner, int *fixedValue, const t_procData * const procData, const int wallIdx);

int p_computeCellOffset(const int outer, const int inner, const int fixedValue, int const * const xlength, const int wallIdx);
#endif
