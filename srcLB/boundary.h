#include "LBDefinitions.h"
#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

/* handles the boundaries in our simulation setup */
void treatBoundary(double *collideField, const t_flagField * const flagField,
	const t_boundPara * const boundPara, const int * const xlength);

#endif
