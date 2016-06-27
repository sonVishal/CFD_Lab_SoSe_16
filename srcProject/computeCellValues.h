#ifndef _COMPUTECELLVALUES_H_
#define _COMPUTECELLVALUES_H_
#include "LBDefinitions.h"

/** computes the number density from the particle distribution functions stored at currentCell.
 *  currentCell thus denotes the address of the first particle distribution function of the
 *  respective cell. The result is stored in density.
 */
void c_computeNumDensity(const double *const currentCell, double *numDensity);

/** computes the velocity within currentCell and stores the result in velocity */
void c_computeVelocity(const double *const currentCell, const double *const c_density,double *c_velocity, const double *c_mass);

/*computes the total velociy as if there were no interacting forces*/
void computeCommonVelocity(const double *const c_density, const double * const c_velocity, const double *const c_tau, double* velocityNI);

/*computes interacting forces between species*/
void c_computeForces(int currentCellIndex, const t_component *const c, const int * const flagField,
    double G[numComp], int * xlength, double forces[3]);

// Computes the equilibrium velocity for each component
void c_computeEqVelocity(const t_component * const c, const double * const commonVelocity,
    const double compDenstiy, const double * const compForce, double compEqVelocity[3]);

/** computes the equilibrium distributions for all particle distribution functions of one
 *  cell from density and velocity and stores the results in feq.
 */
void c_computeFeq(const double density, const double * const velocity, double *feq);

#endif
