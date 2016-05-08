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

    int x,y,z;
    for ( z = 0; z <= xlength+1; ++z) {
        for ( y = 0; y <= xlength+1; ++y) {
            for ( x = 0; x <= xlength+1; ++x) {
                for (int i = 0; i < Q; ++i) {
                    collideField[i] = LATTICEWEIGHTS[i];
                    streamField[i]  = LATTICEWEIGHTS[i];
                }
            }
        }
    }


    /*Lopping over boundary of flagFields*/
    //All points set to zero at initialisation

    int flag  = 1;
    int flagz = 1;

    // b: boundary variable
    for (int b = 0; b <= xlength+1; b = b+xlength+1) {

        if(b > 0){
            flagz = 2;
        }

        for (int i = 0; i < xlength+1; ++i) {
            for (int j = 0; j < xlength+1; ++j) {
                // Remember: 
                //flagField(x,y,z)=flagField[x + y*xlength + z*xlength*xlength]

                /*Set boundary flag at z-boundaries*/
                flagField[j + i*xlength + b*xlength*xlength] = flagz;

                /*Set boundary flag at y-boundaries*/
                flagField[j + b*xlength + i*xlength*xlength] = flag;

                /*Set boundary flag at x-boundaries*/
                flagField[b + j*xlength + i*xlength*xlength] = flag;

            }
        }
    }

}
