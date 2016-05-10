#ifndef _MAIN_C_
#define _MAIN_C_

#include "boundary.h"
#include "collision.h"
#include "initLB.h"
#include "streaming.h"
#include "visualLB.h"
#include "LBDefinitions.h"

int main(int argc, char *argv[]){

    double *collideField=NULL;
    double *streamField=NULL;

    int *flagField=NULL;
    int xlength;
    double tau;
    double velocityWall[3];
    int t = 0;
    int timesteps;
    int timestepsPerPlotting;

    // File printing parameters
    char fName[80];
    snprintf(fName, 80, "pv_files/test");

    /*Read parameters and check the bounds on tau*/
    //tau is calculated automatically from the reynoldsnumber
    int err_check = readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting,argc, &argv[1]);
    if(err_check == -1){
        printf("ERROR: tau<0.5\n");
        return(-1);
    }
    else if(err_check == -2){
        printf("ERROR: tau>2\n");
        return(-2);

    }

    /*Initializing pointers*/
    size_t totalsize = (xlength+2)*(xlength+2)*(xlength+2);
    collideField  = (double *)  malloc(Q*totalsize * sizeof( double ));
    streamField   = (double *)  malloc(Q*totalsize * sizeof( double ));

    /* calloc: only required to set boundary values. Sets every value to zero*/
    flagField     = (int *)  calloc(totalsize, sizeof( int ));

    initialiseFields(collideField, streamField, flagField, xlength);

    printf("INFO: write vtk file at time t = %d \n", t);
    writeVtkOutput(collideField,flagField,fName,t,xlength);

    for(t = 1; t <= timesteps; t++){
	    double *swap=NULL;
	    doStreaming(collideField,streamField,flagField,xlength);

	    swap = collideField;
	    collideField = streamField;
	    streamField = swap;

	    doCollision(collideField,flagField,&tau,xlength);
	    treatBoundary(collideField,flagField,velocityWall,xlength);

	    if (t%timestepsPerPlotting==0){
            printf("INFO: write vtk file at time t = %d \n", t);
	        writeVtkOutput(collideField,flagField,fName,t,xlength);
	    }
    }

    free(streamField);
    free(collideField);
    free(flagField);

    return 0;
}

#endif
