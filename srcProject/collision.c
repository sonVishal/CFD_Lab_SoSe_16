#include "collision.h"

// Get the post collision cell distribution using BGK update rule
// for the current cell
void c_computePostCollisionDistributions(double *currentCell, const double tau, const double *const feq){
	//double tau_val = *tau;
	for ( int i = 0;  i< Q; ++i ) {
		currentCell[i] = currentCell[i]  - (currentCell[i]  - feq[i])/(tau);

#ifndef NDEBUG
		if(currentCell[i]<0){
			char msg[100];
			sprintf(msg, "A negative cell particle distribution (value=%f) was detected!! (Aborting)", currentCell[i]);
			ERROR(msg);
		}
#endif
	}
}

//TODO: (TKS) Adapt to multiple components
// Perform collision for all inner cells
void doCollision(t_component *c, double G[numComp][numComp], int *xlength){

	// Define iteration indices
	int cellidx, idx, x, y, z, n;

	double *currentCell;

	// Perform collision on all "inner" (FLUID) cells
	for (z = 1; z <= xlength[2] ; z++) { //z
		for (y = 1; y <= xlength[1]; y++) { //y
			for (x = 1; x <= xlength[0]; x++) { //x

                // Allocate memory to local cell parameters
                //TODO: (TKS) remove tau array and input c into function instead.
                double c_tau[numComp];              //Array which holds each tau value.
                double c_density[numComp];          //Array with the densities for each component in current cell

                //TODO: (TKS) need to save number density in every cell to use it in the force function.
                double c_numDensity[numComp];       //Array with the number densities for each component
                double c_velocity[numComp];         //Component velocity
                double c_velocityEq[3];    //Array with the velocities for each component
				double c_force[3];

                double commonVelocity[3];  // Total equilibrium velocity if there were no forces between species
                double feq[19];        // Equilibrium distribution of number density

				cellidx = p_computeCellOffsetXYZ(x, y, z, xlength);
				idx = Q*cellidx;

                for (n = 0; n < numComp; ++n){ //c //TODO: (TKS) Possible to optimize this loop nest. We now lose locality (?)

                    // Get the index of the first distribution in current component in the current cell
                    currentCell = &c[n].collideField[idx];

                    // Compute the cell density and number density for each component
                    //TODO: (TKS) need to compute number density in each neighbouring cell.
                    c_computeNumDensity(currentCell, &c_numDensity[n]);
                    c_density[n] = c_numDensity[n]*c[n].m;

#ifndef NDEBUG
                    // We check if the density deviation is more than densityTol%
                    // This value can be changed in LBDefinitions.h
                    if(fabs(c_density[n]-1) > densityTol){
                        char msg[120];
                        sprintf(msg, "A density value (%f) outside the given tolerance of %.2f %% was detected in cell: "
                                "x=%i, y=%i, z=%i, in component %d", c_density[n], (densityTol*100), x, y, z, n);
                        ERROR(msg);
                    }
#endif
                    // Compute the cell velocity for current component
                    c_computeVelocity(currentCell, &c_density[n], &c_velocity[n], &c[n].m);
                }

                computeCommonVelocity(c_density, c_velocity, c_tau, commonVelocity);

				for (n = 0; n < numComp; n++) {
					//Compute the force.
					c_computeForces(cellidx, &c[n], G[n], xlength, c_force);

					//Compute the equilibrium velocity.
					c_computeEqVelocity(&c[n], commonVelocity, c_density[n], c_force, c_velocityEq);

					// Compute the equilibrium distributions
					c_computeFeq(c_density[n], c_velocityEq, feq);

					// Compute the post collision distributions
					c_computePostCollisionDistributions(&c[n].collideField[idx], c[n].tau, feq);
				}
			}
		}
	}
}
