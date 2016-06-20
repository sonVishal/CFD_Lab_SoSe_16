#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "helper.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>

/** handles the boundaries in our simulation setup */
void treatBoundary(t_component *c, int const * const flagField, const t_procData * const procData);

#endif
