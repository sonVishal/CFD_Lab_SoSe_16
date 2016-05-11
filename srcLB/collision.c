#include "collision.h"
#include "LBDefinitions.h"
#include "helper.h"

#include <stdio.h>

// Get the post collision cell distribution using BGK update rule
// for the current cell
void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
    //double tau_val = *tau;
    for ( int i = 0;  i< Q; ++i ) {
        currentCell[i]  = currentCell[i]  - (currentCell[i]  - feq[i])/(*tau);
        //TODO: (TKS) This makes it run 6 sec slower on O3
        if(currentCell[i]<0){
            ERROR("ERROR: Encountered negative particle distribution (Aborting)\n");
        }
    }
    //TODO: (TKS) loops run better for O3, rolling out better for O0
    //      Need to decide which version we are going for.
    //      If going for the one below, need to implement neg- particle check for it.
        //currentCell[0]  = currentCell[0]  - (currentCell[0]  - feq[0])/(*tau);
        //currentCell[1]  = currentCell[1]  - (currentCell[1]  - feq[1])/(*tau);
        //currentCell[2]  = currentCell[2]  - (currentCell[2]  - feq[2])/(*tau);
        //currentCell[3]  = currentCell[3]  - (currentCell[3]  - feq[3])/(*tau);
        //currentCell[4]  = currentCell[4]  - (currentCell[4]  - feq[4])/(*tau);
        //currentCell[5]  = currentCell[5]  - (currentCell[5]  - feq[5])/(*tau);
        //currentCell[6]  = currentCell[6]  - (currentCell[6]  - feq[6])/(*tau);
        //currentCell[7]  = currentCell[7]  - (currentCell[7]  - feq[7])/(*tau);
        //currentCell[8]  = currentCell[8]  - (currentCell[8]  - feq[8])/(*tau);
        //currentCell[9]  = currentCell[9]  - (currentCell[9]  - feq[9])/(*tau);
        //currentCell[10] = currentCell[10] - (currentCell[10] - feq[10])/(*tau);
        //currentCell[11] = currentCell[11] - (currentCell[11] - feq[11])/(*tau);
        //currentCell[12] = currentCell[12] - (currentCell[12] - feq[12])/(*tau);
        //currentCell[13] = currentCell[13] - (currentCell[13] - feq[13])/(*tau);
        //currentCell[14] = currentCell[14] - (currentCell[14] - feq[14])/(*tau);
        //currentCell[15] = currentCell[15] - (currentCell[15] - feq[15])/(*tau);
        //currentCell[16] = currentCell[16] - (currentCell[16] - feq[16])/(*tau);
        //currentCell[17] = currentCell[17] - (currentCell[17] - feq[17])/(*tau);
        //currentCell[18] = currentCell[18] - (currentCell[18] - feq[18])/(*tau);
}

// Perform collision for all inner cells
void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){

    // Define iteration indices
    int idx, x, y, z;

    // Temporary variables for xlength^2 
    long int const xlen2 = (xlength+2)*(xlength+2);

    // Temporary variables for z and y offsets
    int zOffset, yOffset;

    // Perform collision on all "inner" cells
    for (z = 1; z <= xlength ; z++) {
        zOffset = z*xlen2;
        for (y = 1; y <= xlength; y++) {
            yOffset = y*(xlength+2);
            for (x = 1; x <= xlength; x++) {

                // Get the index of the first distribution
                // in the current cell
                idx = Q*(zOffset + yOffset + x);
                double *currentCell = &collideField[idx];

                // Allocate memory to local cell parameters
                double density;
                double velocity[3] = {0,0,0};
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
