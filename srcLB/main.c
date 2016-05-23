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
    double *collideField    =NULL;
    double *streamField     =NULL;
    int *flagField          =NULL;

    // Simulation parameters
    int xlength[3];
    double tau;
    double velocityWall[3];
    int t = 0;
    int timesteps;
    int timestepsPerPlotting;
    t_boundPara boundPara[NUM_WALLS];

    //To be safe allocating memory for max line length (defined in helper.h)
    char problem[MAX_LINE_LENGTH];

    //Timing variables:
    clock_t begin_timing, end_timing;
    double time_spent;

    /*Read parameters and check the bounds on tau*/
    //tau is calculated automatically from the reynoldsnumber
    readParameters(xlength, &tau, boundPara, &timesteps, &timestepsPerPlotting,
    		problem, argc, &argv[1]);

    /* TODO: (DL) the current name includes also the file ending '.pgm' maybe we
     * dont want that..
     */
    // File printing parameters
    char fName[MAX_LINE_LENGTH+9]; // 9 chars for 'pv_files/'
    snprintf(fName, MAX_LINE_LENGTH+9, "pv_files/%s", problem);

#ifdef NO_CHECKS
    printf("INFO: The compiler directive NO_CHECKS is enabled. Faster execution time is gained, "
    		"at the cost of less correctness checks during runtime!\n");
#else
    printf("INFO: The compiler directive NO_CHECKS is disabled. Checks for "
    		"correctness are carried out at the cost of execution speed!\n"
            "      Use \"make speed\" for a faster execution time.\n");
#endif

    /*Allocate memory to pointers*/
    //                   X             Y              Z
    int totalsize = (xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2);

    collideField  = (double *)  malloc(Q*totalsize * sizeof( double ));
    if (collideField == 0) ERROR("Storage cannot be allocated");

    streamField   = (double *)  malloc(Q*totalsize * sizeof( double ));
    if (streamField == 0) ERROR("Storage cannot be allocated");

    /* calloc: only required to set boundary values. Sets every value to zero*/
    flagField     = (int *)  calloc(totalsize, sizeof( int ));
    if (flagField == 0) ERROR("Storage cannot be allocated");

    // Initialize all the fields
    initialiseFields(collideField, streamField, flagField, xlength, boundPara, problem);

    //TODO: (DL) DELETE RETURN VALUE WHEN FINISHED THE INITIALIZATION FOR WS3

    printf("\nINFO: Storing cell data in VTK files.\n      Please use the"
    " \"Cell Data to Point Data\" filter in paraview to view nicely interpolated data. \n\n");

    // Write the VTK at t = 0
    printf("INFO: write vtk file at time t = %d \n", t);
    writeVtkOutput(collideField,flagField,fName,t,xlength);

    return 1;

    begin_timing = clock();
    for(t = 1; t <= timesteps; t++){
	    double *swap = NULL;
	    doStreaming(collideField,streamField,flagField,xlength);

	    swap = collideField;
	    collideField = streamField;
	    streamField = swap;

	    doCollision(collideField,flagField,&tau, xlength);
	    treatBoundary(collideField,flagField,velocityWall, xlength);

	    if (t%timestepsPerPlotting == 0){
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
    printf("Mega Lattice Updates per Seconds: \t %f MLUPS \n",
    		(totalsize/(1000000*time_spent))*timesteps);

    free(streamField);
    free(collideField);
    free(flagField);

    return 0;
}

#endif
