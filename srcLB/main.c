#ifndef _MAIN_C_
#define _MAIN_C_

#include "boundary.h"
#include "collision.h"
#include "initLB.h"
#include "streaming.h"
#include "visualLB.h"
#include "LBDefinitions.h"
#include <time.h>
#include "unitTest.h"


//TODO: (TKS) temporarily have the streamingcollide func here.
//      * Should rename stream and collidefield

void streamCollide(int xlength, double* streamField, double* collideField, double* feq, int* flagField){

	// Define iteration indices
	int cellIdx, fieldIdx;

	// Perform collision on all "inner" (FLUID) cells
	for (int z = 1; z <= xlength ; z++) {
		for (int y = 1; y <= xlength; y++) {
			for (int x = 1; x <= xlength; x++) {

				// Get the index of the first distribution
				// in the current cell
				fieldIdx = p_computeCellOffsetXYZ(x, y, z, xlength);
				cellIdx = Q*fieldIdx;

                
                for (int i = 0; i < Q; ++i) {
                    
                //TODO: Do they stream to other cells in ref?
                //          * We want to take distributions in.
                collideField[cellIdx] = streamField[cellIdx];


                }
            }
        }
    }




}

int main(int argc, char *argv[]){

    // Distribution function vectors
    t_component c[NUMCOMP];

    c[0].tau = 1.0;
    c[0].m = 1.0;
    c[0].psiFctCode = 0;


    // double G[2][2] = {{0.01,0.04},{0.04,0.02}};
    //double G[1][1] = {{-0.27}};

    int *flagField = NULL;
    double *feq = NULL;

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
    for (int i = 0; i < NUMCOMP; i++) {
        c[i].collideField = (double *)  malloc(Q*totalsize * sizeof( double ));
        c[i].streamField = (double *)  malloc(Q*totalsize * sizeof( double ));
        c[i].rho = (double *)  malloc(totalsize * sizeof( double ));
    }

    //TODO: Beware! feq has changed to be of size Q*totalsize.
    feq = (double *)  malloc(Q*totalsize * sizeof( double ));

    initializeUnitTest(totalsize);

    /* calloc: only required to set boundary values. Sets every value to zero*/
    flagField     = (int *)  calloc(totalsize, sizeof( int ));

    // Initialize all the fields
    printf("init fields\n");
    initialiseFields(c, flagField,feq, xlength);
    printf("init fields done\n");

    //return 0;

    printf("\nINFO: Storing cell data in VTK files.\n      Please use the"
    " \"Cell Data to Point Data\" filter in paraview to view nicely interpolated data. \n\n");

    // Write the VTK at t = 0
    printf("INFO: write vtk file at time t = %d \n", t);
    writeVtkOutput(c,flagField,fName,t,xlength);
    // writeVtkOutputDebug(c,flagField,fName,t,xlength);

    tau = 1.0;

    begin_timing = clock();



    //TODO: Wih reference to the ref solution the steps needed to restructure are:
    // Also think of that ref did not use ghost layers.
    //
    // Keeping the names for now to not change notation for now.
    // f     --> stream
    // f_new --> collide
    // swap  --> moved to the end.
    // create combined stream/collide function
    // calc_dPdt             --> computeForce
    // updateDensityVelocity --> computeDensity + computeVelocity
    // updateFeq             --> computeFeq
    
    // treatBoundary inserted where the BC is treated.



    for(t = 1; t <= timesteps; t++){

        treatBoundary(c,flagField,velocityWall,xlength, 1);

        //TODO: Finish this func.
        streamCollide(xlength, c[0].streamField, c[0].collideField, feq, flagField);

        //TODO: compute force
        
        //TODO: compute density
        //TODO: compute velocity
        //TODO: compute feq
        //TODO: streamField = collideField




	    if (t%timestepsPerPlotting == 0){
            printf("INFO: write vtk file at time t = %d \n", t);
	        writeVtkOutput(c,flagField,fName,t,xlength);
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

    for (int i = 0; i < NUMCOMP; i++) {
        free(c[i].streamField);
        free(c[i].collideField);
        free(c[i].rho);
    }
    free(flagField);

    return 0;
}

#endif
