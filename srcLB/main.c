#ifndef _MAIN_C_
#define _MAIN_C_

#include "streamCollide.h"
#include "boundary.h"
#include "initLB.h"
#include "visualLB.h"
#include "LBDefinitions.h"
#include <time.h>
#include "unitTest.h"

//TODO: (TKS) We assume that mass is 1 and use numDensity as density in different places.
//            This must be fixed if we are doing several components.
//TODO: (TKS) Rename streamField and collideField.
//TODO: (TKS) Read in parameters for components from file.
//TODO: (TKS) Comment the code.

int main(int argc, char *argv[]){

    // Distribution function vectors
    t_component c[NUMCOMP];

    c[0].tau = 1.01;
    c[0].m = 1.0;
    c[0].psiFctCode = 0;

    // double G[2][2] = {{0.01,0.04},{0.04,0.02}};
    double G[1][1] = {{-0.27}};

    int *flagField = NULL;

    // Simulation parameters
    int xlength;
    double tau;
    int t = 0;
    int timesteps;
    int timestepsPerPlotting;

    // File printing parameters
    char fName[80];
    snprintf(fName, 80, "pv_files/project");

    //Timing variables:
    clock_t begin_timing, end_timing;
    double time_spent;

    /*Read parameters and check the bounds on tau*/
    readParameters(&xlength, &tau, &timesteps, &timestepsPerPlotting,argc, &argv[1]);

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
    for (int i = 0; i < NUMCOMP; i++) {
        c[i].collideField = (double *)  malloc(Q*totalsize * sizeof( double ));
        c[i].streamField  = (double *)  malloc(Q*totalsize * sizeof( double ));
        c[i].feq          = (double *)  malloc(Q*totalsize * sizeof( double ));
        c[i].rho          = (double *)  malloc(totalsize * sizeof( double ));
        c[i].velocity     = (double **)  malloc(totalsize * sizeof( double* ));
        c[i].force        = (double **)  malloc(totalsize * sizeof( double* ));

        for(int j = 0; j < totalsize; ++j){
            c[i].velocity[j] = (double*) malloc(3* sizeof( double));
            c[i].force[j]    = (double*) malloc(3* sizeof( double));
        }
    }


    initializeUnitTest(totalsize);

    /* calloc: only required to set boundary values. Sets every value to zero*/
    flagField     = (int *)  calloc(totalsize, sizeof( int ));

    // Initialize all the fields
    initialiseFields(c, flagField, xlength);

    printf("\nINFO: Storing cell data in VTK files.\n      Please use the"
    " \"Cell Data to Point Data\" filter in paraview to view nicely interpolated data. \n\n");

    // Write the VTK at t = 0
    printf("INFO: write vtk file at time t = %d \n", t);
    writeVtkOutput(c,flagField,fName,t,xlength);
    // writeVtkOutputDebug(c,flagField,fName,t,xlength);

    begin_timing = clock();


    // Current mapping
    // f     --> stream
    // f_new --> collide
    for(t = 1; t <= timesteps; t++){

        treatBoundary(c,xlength);

        streamCollide(c, xlength, flagField);

		computeForce(c, xlength, flagField, G);

        computeDensityAndVelocity(c, xlength);

        computeFeq(c, &xlength);


        for (int k = 0; k < Q*totalsize; ++k) {
            c[0].streamField[k] = c[0].collideField[k];
        }

        //TODO: (TKS) Swap does not work! Possible to find a way that it does?
        //double *swap = NULL;
        //for (int k = 0; k < NUMCOMP; ++k) {
            //swap              = c[k].streamField;
            //c[k].streamField  = c[k].collideField;
            //c[k].collideField = swap;
        //}


	    if (t%timestepsPerPlotting == 0){
            printf("INFO: write vtk file at time t = %d \n", t);
	        writeVtkOutput(c,flagField,fName,t,xlength);
            timestepsPerPlotting = 2*timestepsPerPlotting;
            // writeVtkOutputDebug(c,flagField,fName,t,xlength);
	    }
        printf("Time t = %d done\n",t);
    }



    end_timing = clock();
    time_spent = (double)(end_timing - begin_timing) / CLOCKS_PER_SEC;

    printf("\n===============================================================\n");
    printf("\nINFO TIMING:\n");
    printf("Execution time (main loop): \t\t %.3f seconds \n", time_spent);
    printf("#cells (including boundary): \t\t %i cells \n", totalsize);
    printf("Mega Lattice Updates per Seconds: \t %f MLUPS \n",
    		(totalsize/(1000000*time_spent))*timesteps);

    //TODO: Remember to free memory
    for (int i = 0; i < NUMCOMP; i++) {
        free(c[i].streamField);
        free(c[i].collideField);
        free(c[i].rho);
        free(c[i].feq);

        for(int j = 0; j < totalsize; ++j){
            c[i].velocity[j] = (double*) malloc(3* sizeof( double));
            c[i].force[j]    = (double*) malloc(3* sizeof( double));
        }
        free(c[i].velocity);
        free(c[i].force);
    }
    free(flagField);

    return 0;
}

#endif
