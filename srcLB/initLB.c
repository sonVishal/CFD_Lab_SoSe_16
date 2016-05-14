
#include <math.h>
#include "initLB.h"
#include "LBDefinitions.h"
#include "helper.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){

    double xvelocity, yvelocity, zvelocity;
    double Re;
    double u_wall;

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
    u_wall = sqrt(xvelocity*xvelocity + yvelocity*yvelocity+zvelocity*zvelocity); 
    *tau    =  u_wall*(*xlength)/(C_S*C_S*Re) +0.5;
    printf("INFO: Calculated tau = %f\n", *tau);

    if(*tau<0.5 || *tau>2){
        ERROR("tau is out of stability region (aborting) \n");
    }

    printf("INFO: Wall speed = %f \n\n", u_wall);

    if(u_wall >= C_S){
    	ERROR("Wall speed is supersonic (aborting). \n");
    }

  return 0;
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){

    /*Setting initial distributions*/
    //f_i(x,0) = f^eq(1,0,0) = w_i

    // current cell index
    int idx;

    // Temporary variables for xlength^2
    int const xlen2 = (xlength+2)*(xlength+2);

    // Temporary variables for z and y offsets
    int zOffset, yzOffset;

    int x,y,z;
    for ( z = 0; z <= xlength+1; ++z) {
        zOffset = z*xlen2;
        for ( y = 0; y <= xlength+1; ++y) {
            yzOffset = y*(xlength+2) + zOffset;
            for ( x = 0; x <= xlength+1; ++x) {
                // Compute the base index
                idx = Q*(yzOffset + x);
                for (int i = 0; i < Q; ++i) {
                    collideField[idx+i] = LATTICEWEIGHTS[i];
                    streamField[idx+i]  = LATTICEWEIGHTS[i];
                }
            }
        }
    }


    /*Lopping over boundary of flagFields*/
    //All points set to zero at initialisation

    // These are the no-slip walls
    // z = 0
    for (y = 0; y <= xlength+1; y++) {
        idx = y*(xlength+2);
        for (x = 0; x <= xlength+1; x++) {
            flagField[x+idx] = 1;
        }
    }

    // x = 0
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2;
        for (y = 0; y <= xlength+1; y++) {
            flagField[zOffset+y*(xlength+2)] = 1;
        }
    }

    // x = xlength+1
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2 + xlength + 1;
        for (y = 0; y <= xlength+1; y++) {
            flagField[zOffset+y*(xlength+2)] = 1;
        }
    }

    // y = 0
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2;
        for (x = 1; x <= xlength; x++) {
            flagField[zOffset+x] = 1;
        }
    }
    // y = xlength+1
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2 + (xlength+1)*(xlength+2);
        for (x = 1; x <= xlength; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    // This is the moving wall
    // z = xlength+1
    zOffset = (xlength+1)*xlen2;
    for (y = 0; y <= xlength+1; y++) {
        idx = zOffset + y*(xlength+2);
        for (x = 0; x <= xlength+1; x++) {
            flagField[x+idx] = 2;
        }
    }
}
