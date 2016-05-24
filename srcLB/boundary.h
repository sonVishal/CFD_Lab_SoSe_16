#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "LBDefinitions.h"

/* handles the boundaries in our simulation setup */
void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity,int *xlength);

/*sets inflow condition*/
void p_handleInflow(int x, int y, int z, int *xlength, t_boundPara *boundPara,
                    double *collideField, const int currentCellIndex);

#endif

