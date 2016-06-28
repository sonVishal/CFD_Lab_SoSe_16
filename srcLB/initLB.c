#include "initLB.h"

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
    int const xlen2 = (xlength+2)*(xlength+2);

    // Temporary variables for z and y offsets
    int zOffset, yzOffset;

    /* initialize collideField and streamField */
    int x,y,z;

    for (int k = 0; k < 2; k++) {
        for ( z = 0; z <= xlength+1; ++z) {
            zOffset = z*xlen2;
            for ( y = 0; y <= xlength+1; ++y) {
                yzOffset = y*(xlength+2) + zOffset;
                for ( x = 0; x <= xlength+1; ++x) {
                    // Compute the base index
                    idx = Q*(yzOffset + x);
                    for (int i = 0; i < Q; ++i) {
                        c[k].collideField[idx+i] = LATTICEWEIGHTS[i];
                        c[k].streamField[idx+i]  = LATTICEWEIGHTS[i];
                    }
                }
            }
        }
    }

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
        zOffset = z*xlen2;
        for (y = 0; y <= xlength+1; y++) {
            flagField[zOffset+y*(xlength+2)] = 1;
        }
    }

    //fixed: x = xlength+1
    //We start at 1 to not include previous cells again from z = 0
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2 + xlength + 1;
        for (y = 0; y <= xlength+1; y++) {
            flagField[zOffset+y*(xlength+2)] = 1;
        }
    }

    //fixed: y = 0
    //from 1:xlength only, to not include cells at upper, lower, left and right edges
    //The edge cells are set in the other loops.
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2;
        for (x = 1; x <= xlength; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    //fixed: y = xlength+1
    //same reasoning for index range as in fixed y=0
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2 + (xlength+1)*(xlength+2);
        for (x = 1; x <= xlength; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    // This is the moving wall. All cells at z=xlength+1 are included (also the edge cells).
    // fixed: z = xlength+1
    zOffset = (xlength+1)*xlen2;
    for (y = 0; y <= xlength+1; y++) {
        idx = zOffset + y*(xlength+2);
        for (x = 0; x <= xlength+1; x++) {
            flagField[x+idx] = 2;
        }
    }
}
