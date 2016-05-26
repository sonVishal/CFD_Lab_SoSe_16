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
void doCollision(double *collideField, t_flagField *flagField,const double * const tau, int *xlength){

	// Define iteration indices
	int idx, x, y, z;

	// Temporary variables for xlength^2
	int const xylen2 = (xlength[0]+2)*(xlength[1]+2);

	// Temporary variables for z and y offsets
	int zOffset, yOffset;

	// Perform collision on all "inner" (FLUID) cells
	for (z = 1; z <= xlength[2] ; z++) {
		zOffset = z*xylen2;
		for (y = 1; y <= xlength[1]; y++) {
			yOffset = y*(xlength[0]+2);
			for (x = 1; x <= xlength[0]; x++) {

				// Get the index of the first distribution
				// in the current cell
				idx = Q*(zOffset + yOffset + x);
				if (flagField[zOffset + yOffset + x].type == FLUID) {
					double *currentCell = &collideField[idx];

					// Allocate memory to local cell parameters
					double density;
					double velocity[3];
					double feq[19];

					// Compute the cell density
					computeDensity(currentCell, &density);

#ifndef NO_CHECKS
					// We check if the density deviation is more than densityTol%
					// This value can be changed in LBDefinitions.h
					if(fabs(density-1) > densityTol){
						char msg[120];
						sprintf(msg, "A density value (%f) outside the given tolerance of %.2f %% was detected in cell: "
								"x=%i, y=%i, z=%i", density, (densityTol*100), x, y, z);
						ERROR(msg);
					}
#endif

					// Compute the cell velocity
					computeVelocity(currentCell, &density, velocity);

					// Compute the equilibrium distributions
					computeFeq(&density,velocity,feq);

					// Compute the post collision distributions
					computePostCollisionDistributions(currentCell,tau,feq);
				}
			}
		}
	}
}
