#include "LBDefinitions.h"

#ifndef _COMPUTECELLVALUES_H_
#define _COMPUTECELLVALUES_H_

/** computes the density from the particle distribution functions stored at currentCell.
 *  currentCell thus denotes the address of the first particle distribution function of the
 *  respective cell. The result is stored in density.
 */
void computeNumDensity(const double *const currentCell, double *density);

/** computes the velocity within currentCell and stores the result in velocity */
void computeVelocity(const double *const currentCell, const double * const density,double *velocity);

/** computes the equilibrium distributions for all particle distribution functions of one
 *  cell from density and velocity and stores the results in feq.
 */
void computeFeq(const double * const density, const double * const velocity, double *feq);

void computeCommonVelocity(const double *const c_density, double c_velocity[2][3], t_component *c, double* commonVel);

void computeForce(const int x, const int y, const int z, const int currentCompIndex,
    const t_component *const c, const int * const flagField,
    double const*const G, int xlength, double forces[3]);

void computeEqVelocity(const t_component * const c, const double * const commonVelocity, const double compDenstiy, const double * const compForce, double compEqVelocity[3]);

#endif
