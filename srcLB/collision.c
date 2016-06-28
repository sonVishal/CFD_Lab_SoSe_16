#include "collision.h"
#include "helper.h"

// Get the post collision cell distribution using BGK update rule
// for the current cell
void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
	//double tau_val = *tau;
	for ( int i = 0;  i< Q; ++i ) {
		currentCell[i] = currentCell[i]  - (currentCell[i]  - feq[i])/(*tau);

#ifndef NO_CHECKS
		if(currentCell[i]<0){
			char msg[100];
			sprintf(msg, "A negative cell particle distribution (value=%f) was detected!! (Aborting)", currentCell[i]);
			ERROR(msg);
		}
#endif
	}
}

// Perform collision for all inner cells
void doCollision(t_component *c, double G[2][2], int *flagField, int xlength){

	// Define iteration indices
	int idx, x, y, z;

	// Temporary variables for xlength^2
	int const xlen2 = (xlength+2)*(xlength+2);

	// Temporary variables for z and y offsets
	int zOffset, yOffset;

	// Perform collision on all "inner" (FLUID) cells
	for (z = 1; z <= xlength ; z++) {
		zOffset = z*xlen2;
		for (y = 1; y <= xlength; y++) {
			yOffset = y*(xlength+2);
			for (x = 1; x <= xlength; x++) {

				// Get the index of the first distribution
				// in the current cell
				idx = Q*(zOffset + yOffset + x);
				double *currentCell = &c[0].collideField[idx];

				// Allocate memory to local cell parameters
				double numDensity[2];
				double density[2];
				double velocity[2][3];
				double commonVel[3];
				double eqVel[3];
				double force[3];
				double feq[19];
				for (int k = 0; k < 2; k++) {
					// Compute the cell numDensity
					computeNumDensity(currentCell, &numDensity[k]);
					density[k] = c[k].m*numDensity[k];
					#ifndef NO_CHECKS
					// We check if the numDensity deviation is more than densityTol%
					// This value can be changed in LBDefinitions.h
					if(fabs(density[k]-1) > densityTol){
						char msg[120];
						sprintf(msg, "A density value (%f) outside the given tolerance of %.2f %% was detected in cell: "
						"x=%i, y=%i, z=%i", density[k], (densityTol*100), x, y, z);
						ERROR(msg);
					}
					#endif
					// Compute the cell velocity
					computeVelocity(currentCell, &numDensity[k], velocity[k]);
				}

				computeCommonVelocity(density, velocity, c, commonVel);

				for (int k = 0; k < 2; k++) {
					computeForce(zOffset + yOffset + x, k, c, flagField, G[k], xlength, force);
					// force[0] = 0.0; force[1] = 0.0; force[2] = 0.0;
					computeEqVelocity(&c[k], commonVel, density[k], force, eqVel);

					computeFeq(&density[k],eqVel,feq);

					computePostCollisionDistributions(currentCell,&c[k].tau,feq);
				}
			}
		}
	}
}
