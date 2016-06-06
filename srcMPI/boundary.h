#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "LBDefinitions.h"

/** handles the boundaries in our simulation setup */
void treatBoundary(double * collideField, int const * const flagField, const t_procData procData);

#endif

