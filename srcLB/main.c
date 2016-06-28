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

    // Distribution function vectors
    t_component c[2];

    c[0].tau = 1.0;
    c[0].m = 1.0;
    c[0].psiFctCode = 0;

    c[1].tau = 0.8;
    c[1].m = 2.0;
    c[1].psiFctCode = 1;

    int *flagField = NULL;

    // Simulation parameters
    int xlength;
    double tau;
    double velocityWall[3];
    int t = 0;
    int timesteps;
    int timestepsPerPlotting;

    // File printing parameters
    char fName[80];
    snprintf(fName, 80, "pv_files/worksheet2");

    //Timing variables:
    clock_t begin_timing, end_timing;
    double time_spent;

    /*Read parameters and check the bounds on tau*/
    //tau is calculated automatically from the reynoldsnumber
    readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting,argc, &argv[1]);

#ifdef NO_CHECKS
    printf("INFO: The compiler directive NO_CHECKS is enabled. Faster execution time is gained, "
    		"at the cost of less correctness checks during runtime!\n");
#else
    printf("INFO: The compiler directive NO_CHECKS is disabled. Checks for "
    		"correctness are carried out at the cost of execution speed!\n"
            "      Use \"make speed\" for a faster execution time.\n");
#endif

    /*Allocate memory to pointers*/
    int totalsize = (xlength+2)*(xlength+2)*(xlength+2);
    c[0].collideField = (double *)  malloc(Q*totalsize * sizeof( double ));
    c[0].streamField = (double *)  malloc(Q*totalsize * sizeof( double ));

    c[1].collideField = (double *)  malloc(Q*totalsize * sizeof( double ));
    c[1].streamField = (double *)  malloc(Q*totalsize * sizeof( double ));

    /* calloc: only required to set boundary values. Sets every value to zero*/
    flagField     = (int *)  calloc(totalsize, sizeof( int ));

    // Initialize all the fields
    initialiseFields(c, flagField, xlength);

    printf("\nINFO: Storing cell data in VTK files.\n      Please use the"
    " \"Cell Data to Point Data\" filter in paraview to view nicely interpolated data. \n\n");

    // Write the VTK at t = 0
    printf("INFO: write vtk file at time t = %d \n", t);
    writeVtkOutput(c,flagField,fName,t,xlength);

    begin_timing = clock();
    for(t = 1; t <= timesteps; t++){
	    double *swap = NULL;
	    doStreaming(c,flagField,xlength);

        for (int i = 0; i < 2; i++) {
            swap = c[i].collideField;
            c[i].collideField = c[i].streamField;
            c[i].streamField = swap;
        }

	    doCollision(c[0].collideField,flagField,&c[0].tau,xlength);
        doCollision(c[1].collideField,flagField,&c[1].tau,xlength);
	    treatBoundary(c,flagField,velocityWall,xlength);

	    if (t%timestepsPerPlotting == 0){
            printf("INFO: write vtk file at time t = %d \n", t);
	        writeVtkOutput(c,flagField,fName,t,xlength);
	    }
    }
    end_timing = clock();
    time_spent = (double)(end_timing - begin_timing) / CLOCKS_PER_SEC;

    printf("\n===============================================================\n");
    printf("\nINFO TIMING:\n");
    printf("Execution time (main loop): \t\t %.3f seconds \n", time_spent);
    printf("#cells (including boundary): \t\t %i cells \n", totalsize);
    printf("Mega Lattice Updates per Seconds: \t %f MLUPS \n",
    		(totalsize/(1000000*time_spent))*timesteps);

    for (int i = 0; i < 2; i++) {
        free(c[i].streamField);
        free(c[i].collideField);
    }
    free(flagField);

    return 0;
}

#endif
