#include "collision.h"
#include "LBDefinitions.h"

// Get the post collision cell distribution using BGK update rule
// for the current cell
void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
    int i;
    for (i = 0; i < Q; i++) {
        currentCell[i] = currentCell[i] - (currentCell[i]-feq[i])/(*tau);
    }
}

// Perform collision for all inner cells
void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){

    // Define iteration indices
    int idx, x, y, z;

    // Temporary variables for xlength^2 and xlength^3
    long int const xlen2 = xlength*xlength;

    // Perform collision on all "inner" cells
    for (z = 1; z <= xlength ; z++) {
        for (y = 1; y <= xlength; y++) {
            for (x = 1; x <= xlength; x++) {

                // Get the index of the first distribution
                // in the current cell
                idx = Q*(z*xlen2 + y*xlength + x);
                double *currentCell = &collideField[idx];

                // Allocate memory to local cell parameters
                double density;
                double velocity[3];
                double feq[19];

                // Compute the cell density
                computeDensity(currentCell, &density);

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
