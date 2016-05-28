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

    // Distribution function arrays
    double *collideField    =NULL;
    double *streamField     =NULL;
    t_flagField *flagField  =NULL;

    // Simulation parameters
    int xlength[3];
    double tau;
    int t = 0;
    int timesteps;
    int timestepsPerPlotting;
    t_boundPara boundPara[NUM_WALLS];

    //To be safe allocating memory for max line length (defined in helper.h)
    char problem[MAX_LINE_LENGTH];

    //Timing variables:
    clock_t begin_timing, end_timing;
    double time_spent;

    /*Read parameters and check for valid settings*/
    readParameters(xlength, &tau, boundPara, &timesteps, &timestepsPerPlotting,
    		problem, argc, &argv[1]);

    // File printing parameters
    char fName[MAX_LINE_LENGTH+9]; // 9 chars for 'pv_files/'
    snprintf(fName, MAX_LINE_LENGTH+9, "pv_files/%s", problem);

#ifdef NDEBUG
    printf("INFO: The compiler directive NDEBUG is enabled. Faster execution time is gained, "
    		"at the cost of less correctness checks during runtime!\n");
#else
    printf("INFO: The compiler directive NDEBUG is disabled. Checks for "
    		"correctness are carried out at the cost of execution speed!\n"
            "      Use \"make speed\" for a faster execution time.\n");
#endif

    /*Allocate memory to pointers*/
    //                   X             Y              Z
    int totalsize = (xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2);

    collideField  = (double *)  malloc(Q*totalsize * sizeof( double ));
    streamField   = (double *)  malloc(Q*totalsize * sizeof( double ));
    flagField     = (t_flagField *)  malloc(totalsize * sizeof( t_flagField ));

    if(! collideField || ! streamField || ! flagField ){
    	ERROR("Storage cannot be allocated");
    }

    memset(flagField, INVALID, totalsize * sizeof(t_flagField)); //set default value to -1 (invalid) for security

    // Initialize all the fields
    initialiseFields(collideField, streamField, flagField, xlength, boundPara, problem);

    printf("\nINFO: Storing cell data in VTK files.\n      Please use the"
    " \"Cell Data to Point Data\" filter in paraview to view nicely interpolated data. \n\n");

    // Write the VTK at t = 0
    printf("INFO: write vtk file at time t = %d \n", t);
    writeVtkOutput(collideField,flagField,fName,t,xlength);
    writeVtkDebug(collideField,flagField,fName,xlength);

    begin_timing = clock();
    for(t = 1; t <= timesteps; t++){
	    double *swap = NULL;
	    doStreaming(collideField,streamField,flagField,xlength);

	    swap = collideField;
	    collideField = streamField;
	    streamField = swap;

	    doCollision(collideField,flagField,&tau, xlength);
	    treatBoundary(collideField,flagField,boundPara,xlength);

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
