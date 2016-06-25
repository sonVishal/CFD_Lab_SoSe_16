#ifndef _COMPUTECELLVALUES_H_
#define _COMPUTECELLVALUES_H_
#include "LBDefinitions.h"

/** computes the number density from the particle distribution functions stored at currentCell.
 *  currentCell thus denotes the address of the first particle distribution function of the
 *  respective cell. The result is stored in density.
 */
void c_computeNumDensity(const double *const currentCell, double *numDensity);

/** computes the velocity within currentCell and stores the result in velocity */
void c_computeVelocity(const double *const currentCell, const double *density,double *velocity, const double *mass);

/** computes the equilibrium distributions for all particle distribution functions of one
 *  cell from density and velocity and stores the results in feq.
 */
void computeFeq(const double density, const double * const velocity, double *feq);

#endif
