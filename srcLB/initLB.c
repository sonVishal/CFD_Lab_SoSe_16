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

    //TODO: (TKS) Do a smarter thing to get the correct ws at the right places.
    //TODO: (TKS) NOT finished
    //double w[3] = {w1, w2, w3};

    int index;
    for (int x = 0; x < xlength+1; ++x) {
        for (int y = 0; y < xlength+1; ++y) {
            for (int z = 0; z < xlength+1; ++z) {
                for (int i = 0; i < Q; ++i) {
                    index = Q*( (z+y)*xlength*xlength + x );
                    printf("index = %d\n", index);
                    //collideField[index] = w1;
                    
                }
            }
        }
    }
    
    printf("Initialize fields\n");
}
