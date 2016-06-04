#ifndef _MAIN_C_
#define _MAIN_C_

#include "boundary.h"
#include "collision.h"
#include "initLB.h"
#include "streaming.h"
#include "visualLB.h"
#include "LBDefinitions.h"
#include <time.h>
#include "debug.h"

int main(int argc, char *argv[]){

	// Distribution function vectors
    double *collideField    =NULL;
    double *streamField     =NULL;
    int *flagField          =NULL;

    // Simulation parameters
    int xlength;
    double tau;
    double velocityWall[3];
    int t = 0;
    int timesteps;
    int timestepsPerPlotting;

    // MPI parameters
    int rank;
    int number_of_ranks;
    int iProc, jProc, kProc;

    // Send and read buffers for all possible directions:
    // Look at enum for index and direction correlation
    double *sendBuffer[6];
    double *readBuffer[6];

    // File printing parameters
    char fName[80];
    snprintf(fName, 80, "pv_files/worksheet4");

    // TODO
    initialiseMPI(&rank,&number_of_ranks,argc,argv);

    //Timing variables:
    clock_t begin_timing, end_timing;
    double time_spent;

    /*Read parameters and check the bounds on tau*/
    //tau is calculated automatically from the reynoldsnumber
    readParameters(&xlength, &tau, velocityWall, &iProc, &jProc, &kProc, &timesteps, &timestepsPerPlotting,argc, &argv[1]);

#ifdef NO_CHECKS
    printf("INFO: The compiler directive NO_CHECKS is enabled. Faster execution time is gained, "
    		"at the cost of less correctness checks during runtime!\n");
#else
    printf("INFO: The compiler directive NO_CHECKS is disabled. Checks for "
    		"correctness are carried out at the cost of execution speed!\n"
            "      Use \"make speed\" for a faster execution time.\n");
#endif

    /*Allocate memory to pointers*/
    // TODO: Allocate memory based on domain decomposition
    int totalsize = (xlength+2)*(xlength+2)*(xlength+2);
    collideField  = (double *)  malloc(Q*totalsize * sizeof( double ));
    streamField   = (double *)  malloc(Q*totalsize * sizeof( double ));

    /* calloc: only required to set boundary values. Sets every value to zero*/
    // TODO: allocate based on rank
    flagField     = (int *)  calloc(totalsize, sizeof( int ));

    // TODO: Initialize all the fields based on the rank
    initialiseFields(collideField, streamField, flagField,
        xlength, rank, number_of_ranks);

    // TODO:
    initialiseBuffers(sendBuffer, readBuffer, xlength);

    printf("\nINFO: Storing cell data in VTK files.\n      Please use the"
    " \"Cell Data to Point Data\" filter in paraview to view nicely interpolated data. \n\n");

    // Write the VTK at t = 0
    printf("INFO: write vtk file at time t = %d \n", t);
    writeVtkOutput(collideField,flagField,fName,t,xlength,rank,number_of_ranks);

    begin_timing = clock();
    for(t = 1; t <= timesteps; t++){
	    double *swap = NULL;
        /* TODO:
         do extraction , swap , injection for x (left to right)
         do extraction , swap , injection for x (right to left)
         do extraction , swap , injection for y (forth and back; back and forth)
         do extraction , swap , injection for z (down and up ; up and down)
        */

	    doStreaming(collideField,streamField,flagField,xlength);

	    swap = collideField;
	    collideField = streamField;
	    streamField = swap;

	    doCollision(collideField,flagField,&tau,xlength);
	    treatBoundary(collideField,flagField,velocityWall,xlength);

	    if (t%timestepsPerPlotting == 0){
            printf("INFO: write vtk file at time t = %d \n", t);
            // TODO: (VS) Change the coordinates based on rank
            // TODO: (VS) Check if the time-rank strategy is okay for paraview
	        writeVtkOutput(collideField,flagField,fName,t,xlength,rank,number_of_ranks);
	    }
    }
    // TODO: As of now this is in initLB - I have no idea where to put this function definition
    finaliseMPI();

    // TODO: The times from all the processes need to be added at the end
    end_timing = clock();
    time_spent = (double)(end_timing - begin_timing) / CLOCKS_PER_SEC;

    printf("\n===============================================================\n");
    printf("\nINFO TIMING:\n");
    printf("Execution time (main loop): \t\t %.3f seconds \n", time_spent);
    printf("#cells (including boundary): \t\t %i cells \n", totalsize);
    printf("Mega Lattice Updates per Seconds: \t %f MLUPS \n",
    		(totalsize/(1000000*time_spent))*timesteps);

    /* TODO: (DL) when checking with MPI implementation, the values of collideField
     * have to get collected from the different processes and then checked with
     * the 'reference' solution.
     */

    /* To generate the current reference solution switch to branch:
     * git checkout generate_reference_solution
     */

    //char fileRef[] = {"debug/collideField"};
    // writeCollideFieldDebug(fileRef, collideField, Q*totalsize);
    // checkCollideFieldDebug(fileRef, collideField, Q*totalsize);

    free(streamField);
    free(collideField);
    free(flagField);

    return 0;
}

#endif
