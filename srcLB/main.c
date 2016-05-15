#ifndef _MAIN_C_
#define _MAIN_C_

#include "boundary.h"
#include "collision.h"
#include "initLB.h"
#include "streaming.h"
#include "visualLB.h"
#include "LBDefinitions.h"
#include <time.h>

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
    snprintf(fName, 80, "pv_files/worksheet2");

    /*Read parameters and check the bounds on tau*/
    //tau is calculated automatically from the reynoldsnumber
    readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting,argc, &argv[1]);

#ifdef NO_CHECKS
    printf("INFO: The compiler directive NO_CHECKS is enabled. Faster execution time is gained, "
    		"at the cost of less correctness checks during runtime!\n");
#else
    printf("INFO: The compiler directive NO_CHECKS is disables. Checks for "
    		"correctness are carried out at the cost of execution speed!\n");
#endif

    /*Initializing pointers*/
    int totalsize = (xlength+2)*(xlength+2)*(xlength+2);
    collideField  = (double *)  malloc(Q*totalsize * sizeof( double ));
    streamField   = (double *)  malloc(Q*totalsize * sizeof( double ));

    /* calloc: only required to set boundary values. Sets every value to zero*/
    flagField     = (int *)  calloc(totalsize, sizeof( int ));

    initialiseFields(collideField, streamField, flagField, xlength);

    printf("INFO: write vtk file at time t = %d \n", t);
    writeVtkOutput(collideField,flagField,fName,t,xlength);

    //Timing variables:
    clock_t begin_timing, end_timing;
    double time_spent;

    begin_timing = clock();
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
    end_timing = clock();
    time_spent = (double)(end_timing - begin_timing) / CLOCKS_PER_SEC;

    printf("\n===============================================================\n");
    printf("\nINFO TIMING:\n");
    printf("Execution time (main loop): \t\t %.3f seconds \n", time_spent);
    printf("#cells (including boundary): \t\t %i cells \n", totalsize);
    printf("Mega Lattice Updates per Seconds: \t %f MLUPS \n", (totalsize*timesteps)/(1000000*time_spent));

    free(streamField);
    free(collideField);
    free(flagField);

    return 0;
}

#endif
