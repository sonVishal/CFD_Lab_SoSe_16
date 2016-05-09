#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){

    double xvelocity, yvelocity, zvelocity;

    READ_INT(*argv, *xlength);
    READ_DOUBLE(*argv, *tau);

    READ_DOUBLE(*argv, xvelocity);
    READ_DOUBLE(*argv, yvelocity);
    READ_DOUBLE(*argv, zvelocity);

    velocityWall[0] = xvelocity;
    velocityWall[1] = yvelocity;
    velocityWall[2] = zvelocity;

    READ_INT(*argv, *timesteps);
    READ_INT(*argv, *timestepsPerPlotting);

  return 0;
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){

    /*Setting initial distributions*/
    //f_i(x,0) = f^eq(1,0,0) = w_i

    // current cell index
    int idx;

    // Temporary variables for xlength^2 and xlength^3
    int const xlen2 = xlength*xlength;

    // Temporary variables for z and y offsets
    int zOffset, yOffset;

    int x,y,z;
    for ( z = 0; z <= xlength+1; ++z) {
        zOffset = z*xlen2;
        for ( y = 0; y <= xlength+1; ++y) {
            yOffset = y*xlength;
            for ( x = 0; x <= xlength+1; ++x) {
                // Compute the base index
                idx = Q*(zOffset + yOffset + x);
                for (int i = 0; i < Q; ++i) {
                    collideField[idx+i] = LATTICEWEIGHTS[i];
                    streamField[idx+i]  = LATTICEWEIGHTS[i];
                }
            }
        }
    }


    /*Lopping over boundary of flagFields*/
    //All points set to zero at initialisation

    int flag  = 1;
    int flagz = 1;

    // Temporary offset variables
    int bOffset1, bOffset2, iOffset1, iOffset2;

    // b: boundary variable
    for (int b = 0; b <= xlength+1; b = b+xlength+1) {
        bOffset1 = b*xlength;
        bOffset2 = b*xlen2;

        if(b > 0){
            flagz = 2;
        }

        for (int i = 0; i < xlength+1; ++i) {
            iOffset1 = i*xlength;
            iOffset2 = i*xlen2;
            
            for (int j = 0; j < xlength+1; ++j) {
                // Remember:
                //flagField(x,y,z)=flagField[x + y*xlength + z*xlength*xlength]

                /*Set boundary flag at z-boundaries*/
                flagField[j + iOffset1 + bOffset2] = flagz;

                /*Set boundary flag at y-boundaries*/
                flagField[j + bOffset1 + iOffset2] = flag;

                /*Set boundary flag at x-boundaries*/
                flagField[b + j*xlength + iOffset2] = flag;

            }
        }
    }

}
