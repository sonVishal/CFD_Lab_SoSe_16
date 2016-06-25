#include "collision.h"

// Get the post collision cell distribution using BGK update rule
// for the current cell
void computePostCollisionDistributions(double *currentCell, const double tau, const double *const feq){
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
void doCollision(t_component *c, const int numComp, double** G, int *xlength){

	// Define iteration indices
	int cellidx, idx, x, y, z, n;

	// Temporary variables for xlength^2
	int const xylen = (xlength[0]+2)*(xlength[1]+2);

	// Temporary variables for z and y offsets
	int zOffset, yOffset;

	// Perform collision on all "inner" (FLUID) cells
	for (z = 1; z <= xlength[2] ; z++) { //z
		zOffset = z*xylen;
		for (y = 1; y <= xlength[1]; y++) { //y
			yOffset = y*(xlength[0]+2);
			for (x = 1; x <= xlength[0]; x++) { //x
                
                // Allocate memory to local cell parameters
                //TODO: (TKS) remove tau array and input c into function instead.
                double c_tau[numComp];              //Array which holds each tau value.
                double c_density[numComp];          //Array with the densities for each component

                //TODO: (TKS) need to save number density in every cell to use it in the force function.
                double c_numDensity[numComp];       //Array with the number densities for each component
                double c_velocity[numComp];         //Component velocity
                double c_velocityEq[numComp][3];    //Array with the velocities for each component

                double density;
                double velocityNI[3];  // Total equilibrium velocity if there were no forces between species
                double feq[19];        // Equilibrium distribution of number density

                for (n=0; n < numComp; ++n){ //c //TODO: (TKS) Possible to optimize this loop nest. We now lose locality (?)

                    // Get the index of the first distribution in current component in the current cell
                    idx = Q*(zOffset + yOffset + x);
                    double *currentCell = &c[n].collideField[idx];

                    // Compute the cell density and number density for each component
                    //TODO: (TKS) need to compute number density in each neighbouring cell. 
                    c_computeNumDensity(currentCell, &c_numDensity[n]);
                    c_density[n] = c_numDensity[n]*c->m;

                    //TODO: (TKS) Fix this test to encorporate multiple components
                    //#ifndef NDEBUG
                    //// We check if the density deviation is more than densityTol%
                    //// This value can be changed in LBDefinitions.h
                    //if(fabs(density-1) > densityTol){
                        //char msg[120];
                        //sprintf(msg, "A density value (%f) outside the given tolerance of %.2f %% was detected in cell: "
                                //"x=%i, y=%i, z=%i", density, (densityTol*100), x, y, z);
                        //ERROR(msg);
                    //}
                    //#endif
                
                    // Compute the cell velocity for current component
                    c_computeVelocity(currentCell, &c_density[numComp], &c_velocity[n], &c->m);
                }

                computeVelocityNI(&numComp, c_density, c_velocity, c_tau, velocityNI);
                //TODO: (TKS) Add forces.
                //TODO: (TKS) Add equilibrium velocity.

                //TODO: (TKS) Compute new equilibrium distribution.
                // Compute the equilibrium distributions
                computeFeq(density,velocity,feq);

                // Compute the post collision distributions
                computePostCollisionDistributions(currentCell,c[n]tau,feq);
			}
		}
	}
}
