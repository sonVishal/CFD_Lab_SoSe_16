#include "LBDefinitions.h"
#include "assert.h"

#ifndef _COMPUTECELLVALUES_H_
#define _COMPUTECELLVALUES_H_

/** computes the density from the particle distribution functions stored at currentCell.
 *  currentCell thus denotes the address of the first particle distribution function of the
 *  respective cell. The result is stored in density.
 */
void computeNumDensity(const double *const currentCell, double *density);

/** computes the velocity within currentCell and stores the result in velocity */
void computeVelocity(const double *const currentCell, const double * const density,double *velocity);


void computeDensityAndVelocity(t_component *c, int xlength);

/** computes the equilibrium distributions for all particle distribution functions of one
 *  cell from density and velocity and stores the results in feq.
 */
void computeFeqCell(const double * const density, const double * const velocity, double *feq);
void computeFeq(t_component *c, const int *xlength);

void computeCommonVelocity(const double *const c_density, double c_velocity[2][3], t_component *c, double* commonVel);

void computeForce(t_component *c, int xlength, int *flagField, double G[NUMCOMP][NUMCOMP]);

void computeEqVelocity(const t_component * const c, const double * const commonVelocity, const double compDenstiy, const double * const compForce, double compEqVelocity[3]);

#endif
