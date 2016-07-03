#ifndef _MAIN_C_
#define _MAIN_C_

#include "streamCollide.h"
#include "boundary.h"
#include "collision.h"
#include "initLB.h"
#include "streaming.h"
#include "visualLB.h"
#include "LBDefinitions.h"
#include <time.h>
#include "unitTest.h"

//TODO: Move this into computeCellValues when finished
void computeForce_new(t_component *c, int xlength, int *flagField, double G[NUMCOMP][NUMCOMP]);
void computeDensityAndVelocity(t_component *c, int xlength);

int main(int argc, char *argv[]){

    // Distribution function vectors
    t_component c[NUMCOMP];

    c[0].tau = 1.0;
    c[0].m = 1.0;
    c[0].psiFctCode = 0;

    // double G[2][2] = {{0.01,0.04},{0.04,0.02}};
    double G[1][1] = {{-0.27}};

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

    //return 0;

    printf("\nINFO: Storing cell data in VTK files.\n      Please use the"
    " \"Cell Data to Point Data\" filter in paraview to view nicely interpolated data. \n\n");

    // Write the VTK at t = 0
    printf("INFO: write vtk file at time t = %d \n", t);
    writeVtkOutput(c,flagField,fName,t,xlength);
    // writeVtkOutputDebug(c,flagField,fName,t,xlength);

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

    double *swap = NULL;
    for(t = 1; t <= 1; t++){

        treatBoundary(c,velocityWall,xlength);

        //TODO: need to include
        streamCollide(c, xlength, flagField);
        printf("collideField 10 @ (41,12,13) %.11f\n",c[0].collideField[Q*(41+12*(xlength+2)+13*(xlength+2)*(xlength+2))+10]);


        // int x = 2;
        // int y= 2;
        // int z = 2;
        // int cell = p_computeCellOffsetXYZ_Q(x, y, z, xlength);
        // int dist = 2;

        //TODO: compute force
		computeForce_new(c, xlength, flagField, G);
        printf("Force x @ (6,6,6) %.16f\n",c[0].force[(6+6*(xlength+2)+6*(xlength+2)*(xlength+2))][0]);
        printf("Force y @ (6,6,6) %.16f\n",c[0].force[(6+6*(xlength+2)+6*(xlength+2)*(xlength+2))][1]);
        printf("Force z @ (6,6,6) %.16f\n",c[0].force[(6+6*(xlength+2)+6*(xlength+2)*(xlength+2))][2]);

        //TODO: compute density & compute velocity
        computeDensityAndVelocity(c, xlength);
        printf("Rho @ (6,6,6) %.16f\n",c[0].rho[(6+6*(xlength+2)+6*(xlength+2)*(xlength+2))]);
        printf("Ux @ (6,6,6) %.16f\n",c[0].velocity[(6+6*(xlength+2)+6*(xlength+2)*(xlength+2))][0]);
        printf("Uy @ (6,6,6) %.16f\n",c[0].velocity[(6+6*(xlength+2)+6*(xlength+2)*(xlength+2))][1]);
        printf("Uz @ (6,6,6) %.16f\n",c[0].velocity[(6+6*(xlength+2)+6*(xlength+2)*(xlength+2))][2]);

        //TODO: compute feq
        updateFeq(c, &xlength);
        printf("F_eq 10 @ (41,12,13) %.11f\n",c[0].feq[Q*(41+12*(xlength+2)+13*(xlength+2)*(xlength+2))+10]);

        //TODO: streamField = collideField
        for (int k = 0; k < totalsize; ++k) {
            swap = c[k].collideField;
            c[k].collideField = c[k].streamField;
            c[k].streamField = swap;
        }

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

    //TODO: Remember to free memory
    for (int i = 0; i < NUMCOMP; i++) {
        free(c[i].streamField);
        free(c[i].collideField);
        free(c[i].rho);
    }
    free(flagField);

    return 0;
}

#endif
