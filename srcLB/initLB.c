#include "initLB.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){

    double xvelocity, yvelocity, zvelocity;
    double Re, u_wall, machNr;

    /* Read values from file given in argv */
    READ_INT(*argv, *xlength);
    READ_DOUBLE(*argv, Re);

    READ_DOUBLE(*argv, xvelocity);
    READ_DOUBLE(*argv, yvelocity);
    READ_DOUBLE(*argv, zvelocity);

    velocityWall[0] = xvelocity;
    velocityWall[1] = yvelocity;
    velocityWall[2] = zvelocity;

    READ_INT(*argv, *timesteps);
    READ_INT(*argv, *timestepsPerPlotting);

    /*Calculates tau from the Reynolds number*/
    u_wall  = sqrt(xvelocity*xvelocity + yvelocity*yvelocity+zvelocity*zvelocity);
    *tau    =  u_wall*(*xlength)/(C_S*C_S*Re)+0.5;
    machNr  = u_wall/C_S;

    printf("\nINFO: Calculated tau = %f \n", *tau);
    printf("\nINFO: Wall speed = %f \n", u_wall);
    printf("\nINFO: Mach number = %f \n\n", machNr);

    /* valid settings check*/
    if(u_wall >= C_S){
    	ERROR("Wall speed is supersonic (aborting)! \n");
    }

    if(*tau<=0.5 || *tau>2){
        ERROR("Tau is out of stability region (0.5,2.0) (aborting)! \n");
    }

    /*We allow user defined mach number tolerance for Ma << 1 (default = 0.1)
      To change please look at LBDefinitions.h*/
    if(machNr > machNrTol){
        char buffer[80];
        snprintf(buffer, 80, "Mach number is larger than %f (aborting)! \n",machNrTol);
    	ERROR(buffer);
    }

  return 0;
}

void initialiseFields(t_component * c, int *flagField, int xlength){

    /*Setting initial distributions*/
    //f_i(x,0) = f^eq(1,0,0) = w_i
    // current cell index
    int idx;

    // Temporary variables for xlength^2
    int const xlen2 = xlength+2;
    int const xlen2sq = xlen2*xlen2;

    // Temporary variables for z and y offsets
    int zOffset, yzOffset;

    /* initialize collideField and streamField */

    //Intermediate values to calculate feq for a random density.
    double v[3] = {0,0,0};          //Assume no initial velocity

    //TODO: (TKS) Maybe set the percentage 0.01 as a variable we can choose instead.
    double rhoVar = 0.01*rhoRef;    //How much initial difference is allowed in density-

    int x,y,z;
    srand(1);
    for (int k = 0; k < NUMCOMP; k++) {
        for ( z = 0; z <= xlength+1; ++z) {
            if (z == xlength/2) {
                srand(3);
            }
            zOffset = z*xlen2sq;
            for ( y = 0; y <= xlength+1; ++y) {
                yzOffset = y*(xlength+2) + zOffset;
                for ( x = 0; x <= xlength+1; ++x) {
                    // Compute the base index
                    idx = Q*(yzOffset + x);
                    double feq[19];
                    //TODO: (TKS) Need to do componentwise if we introduce several components.
                    
                    //Set the initial density to a random offsett to rhoRef
                    double rho = rhoRef - 0.5*rhoVar + rhoVar*((double)rand()/(double)RAND_MAX);
                    computeFeq(&rho, v, feq);
                    for (int i = 0; i < Q; ++i) {
                        c[k].collideField[idx+i] = feq[i];
                        c[k].streamField[idx+i]  = feq[i];
                    }
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
