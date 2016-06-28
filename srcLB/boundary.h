#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "LBDefinitions.h"
#include "helper.h"
#include "computeCellValues.h"
#include <stdio.h>

/** handles the boundaries in our simulation setup */
void treatBoundary(t_component *collideField, int* flagField, const double * const wallVelocity,int xlength);

#endif
