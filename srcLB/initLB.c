#include "initLB.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

int readParameters(int *xlength, double *tau, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){

    /* Read values from file given in argv */
    READ_INT(*argv, *xlength);

    READ_INT(*argv, *timesteps);
    READ_INT(*argv, *timestepsPerPlotting);


    //TODO: (TKS) Adapt this to component case.
    //if(*tau<=0.5 || *tau>2){
        //ERROR("Tau is out of stability region (0.5,2.0) (aborting)! \n");
    //}

  return 0;
}

void initialiseFields(t_component * c, int *flagField, int xlength){

    /*Setting initial distributions*/
    /*Initializes to equilibrium state with a random density with 1% max from reference density*/

    // current cell index
    int idx, cellIdx;

    // Temporary variables for xlength^2
    int const xlen2 = xlength+2;
    int const xlen2sq = xlen2*xlen2;

    // Temporary variables for z and y offsets
    int zOffset, yzOffset;

    /* initialize collideField and streamField */

    //Intermediate values to calculate feq for a random density.

    int x,y,z;
    for ( z = 1; z <= xlength; ++z) {
        zOffset = z*xlen2sq;
        for ( y = 1; y <= xlength; ++y) {
            yzOffset = y*(xlength+2) + zOffset;
            for ( x = 1; x <= xlength; ++x) {
                // Compute the base index
                cellIdx = (yzOffset + x);
                idx = Q*cellIdx;

                //Set the initial density to a random offsett to rhoRef
                double rnd = ((double)rand()/(double)RAND_MAX);
                if (rnd <= 0.6) {
                    c[0].rho[cellIdx] = 1.0;
                    c[1].rho[cellIdx] = 0.0;
                } else {
                    c[0].rho[cellIdx] = 0.0;
                    c[1].rho[cellIdx] = 1.0;
                }

                c[0].velocity[cellIdx][0] = 0.0;c[0].velocity[cellIdx][1] = 0.0;c[0].velocity[cellIdx][2] = 0.0;
                c[1].velocity[cellIdx][0] = 0.0;c[1].velocity[cellIdx][1] = 0.0;c[2].velocity[cellIdx][2] = 0.0;

                computeFeqCell(&c[0].rho[cellIdx], c[0].velocity[cellIdx], &c[0].feq[idx]);
                computeFeqCell(&c[1].rho[cellIdx], c[1].velocity[cellIdx], &c[1].feq[idx]);
                for (int i = 0; i < Q; ++i) {
                    c[0].collideField[idx+i] = c[0].feq[idx +i];
                    c[0].streamField[idx+i]  = c[0].feq[idx +i];
                    c[1].collideField[idx+i] = c[1].feq[idx +i];
                    c[1].streamField[idx+i]  = c[1].feq[idx +i];
                }
            }
        }
    }


    // NOTE: For debug
    // for (z = 1; z <= xlength; z++) {
    //     zOffset = z*xlen2sq;
    //     for (x = 1; x <= xlength; x++) {
    //         idx = Q*(zOffset + 1*(xlength+2) + x);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 1;
    //             c[0].streamField[idx+i] = 1;
    //         }
    //         idx = Q*(zOffset + (xlength)*(xlength+2) + x);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 2;
    //             c[0].streamField[idx+i] = 2;
    //         }
    //     }
    // }
    //
    // for (z = 1; z <= xlength; z++) {
    //     zOffset = z*xlen2sq;
    //     for (y = 1; y <= xlength; y++) {
    //         idx = Q*(zOffset + y*(xlength+2) + 1);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 3;
    //             c[0].streamField[idx+i] = 3;
    //         }
    //         idx = Q*(zOffset + y*(xlength+2) + xlength);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 4;
    //             c[0].streamField[idx+i] = 4;
    //         }
    //     }
    // }
    //
    // for (y = 1; y <= xlength; y++) {
    //     for (x = 1; x <= xlength; x++) {
    //         idx = Q*(xlen2sq + y*(xlength+2) + x);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 5;
    //             c[0].streamField[idx+i] = 5;
    //         }
    //         idx = Q*(xlen2sq*(xlength) + y*(xlength+2) + x);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 6;
    //             c[0].streamField[idx+i] = 6;
    //         }
    //     }
    // }
    // y = 1; z = 1;
    // for (x = 1; x <= xlength; x++) {
    //     idx = Q*(z*xlen2sq + y*xlen2 + x);
    //     for (int i = 0; i < Q; i++) {
    //         c[0].collideField[idx + i] = 1;
    //         c[0].streamField[idx + i] = 1;
    //     }
    // }
    //
    // y = xlength; z = xlength;
    // for (x = 1; x <= xlength; x++) {
    //     idx = Q*(z*xlen2sq + y*xlen2 + x);
    //     for (int i = 0; i < Q; i++) {
    //         c[0].collideField[idx + i] = 3;
    //         c[0].streamField[idx + i] = 3;
    //     }
    // }

    // x= 1; y = 1;
    // for (z = 1; z <= xlength; z++) {
    //     idx = Q*(z*xlen2sq + y*xlen2 + x);
    //     for (int i = 0; i < Q; i++) {
    //         c[0].collideField[idx + i] = 1;
    //         c[0].streamField[idx + i] = 1;
    //     }
    // }
    //
    // x = xlength; y = xlength;
    // for (z = 1; z <= xlength; z++) {
    //     idx = Q*(z*xlen2sq + y*xlen2 + x);
    //     for (int i = 0; i < Q; i++) {
    //         c[0].collideField[idx + i] = 3;
    //         c[0].streamField[idx + i] = 3;
    //     }
    // }

    /*Looping over boundary of flagFields*/
    //All points set to zero at memory allocation (using calloc)

    //These are the no-slip walls
    //fixed: z = 0
    for (y = 0; y <= xlength+1; y++) {
        idx = y*(xlength+2);
        for (x = 0; x <= xlength+1; x++) {
            flagField[x+idx] = 1;
        }
    }

    //fixed: x = 0
    //We start at 1 to not include previous cells again from z = 0
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2sq;
        for (y = 0; y <= xlength+1; y++) {
            flagField[zOffset+y*(xlength+2)] = 1;
        }
    }

    //fixed: x = xlength+1
    //We start at 1 to not include previous cells again from z = 0
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2sq + xlength + 1;
        for (y = 0; y <= xlength+1; y++) {
            flagField[zOffset+y*(xlength+2)] = 1;
        }
    }

    //fixed: y = 0
    //from 1:xlength only, to not include cells at upper, lower, left and right edges
    //The edge cells are set in the other loops.
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2sq;
        for (x = 1; x <= xlength; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    //fixed: y = xlength+1
    //same reasoning for index range as in fixed y=0
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2sq + (xlength+1)*(xlength+2);
        for (x = 1; x <= xlength; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    // This is the moving wall. All cells at z=xlength+1 are included (also the edge cells).
    // fixed: z = xlength+1
    zOffset = (xlength+1)*xlen2sq;
    for (y = 0; y <= xlength+1; y++) {
        idx = zOffset + y*(xlength+2);
        for (x = 0; x <= xlength+1; x++) {
            flagField[x+idx] = 2;
        }
    }
}
