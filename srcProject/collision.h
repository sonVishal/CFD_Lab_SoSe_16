#ifndef _COLLISION_H_
#define _COLLISION_H_
#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "helper.h"

/** computes the post-collision distribution functions according to the BGK update rule and
 *  stores the results again at the same position.
 */
void computePostCollisionDistributions(double *currentCell, const double tau, const double * const feq);


/** carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
void doCollision(t_component *c1, int numComp, double **G, int *xlength);
#endif
