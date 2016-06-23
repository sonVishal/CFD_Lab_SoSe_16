#ifndef _MAIN_C_
#define _MAIN_C_
#include "boundary.h"
#include "collision.h"
#include "initLB.h"
#include "streaming.h"
#include "visualLB.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "parallel.h"
#include "unitTest.h"
#include <stdio.h>
#include <mpi/mpi.h>


//TODO: (TKS) Where should these swap functions be placed?
void swapComponentFields(t_component *c, int numComp){
    for (int i = 0; i < numComp; ++i) {
        double* swap = c[i].collideField;
        c[i].collideField = c[i].streamField;
        c[i].streamField = swap;
    }
}

int main(int argc, char *argv[]){

	// Distribution function vectors
    double *collideField    = NULL;
    double *streamField     = NULL;
    int *flagField          = NULL;

    //TODO: (TKS) Should read in number of components externally. For now having a global variable.
    // Array of the components in our simulation
    t_component c[g_numComp];

    //TODO: Temporarily initialization of component, REMOVE when implementing multicomponent
    c[0].collideField = collideField;
    c[0].streamField  = streamField;

    // Simulation parameters
    int xlength = 0;
    t_procData procData;
    double wallVelocity[3];

    int t = 0;
    int timesteps;
    int timestepsPerPlotting;

    // MPI parameters
    // MPI_Status status;
    int procsPerAxis[3];

    // Send and read buffers for all possible directions:
    // Look at enum for index and direction correlation
    double *sendBuffer[6];
    double *readBuffer[6];

    // initialiseMPI(&procData.rank,&procData.numRanks,argc,argv);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procData.numRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &procData.rank);

    // File printing parameters
    char fName[80];
    snprintf(fName, 80, "pv_files/WS4_rank");

    //Timing variables:
    double beginProcTime, endProcTime, beginSimTime, endSimTime;

    /* Read parameters and check the bounds on tau
     * Only performed by the root and broadcasted*/
    //tau is calculated automatically from the reynoldsnumber
    if (procData.rank == 0) {
        readParameters(&xlength, &c[0].tau, wallVelocity, procsPerAxis, &timesteps, &timestepsPerPlotting, argc, &argv[1]);
    }

    // Broadcast the data from rank 0 (root) to other processes
    MPI_Bcast(&xlength, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&c[0].tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(wallVelocity, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&timestepsPerPlotting, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(procsPerAxis, 3, MPI_INT, 0, MPI_COMM_WORLD);

    // Assign the wall velocity
    procData.wallVelocity[0] = wallVelocity[0];
    procData.wallVelocity[1] = wallVelocity[1];
    procData.wallVelocity[2] = wallVelocity[2];

    // Abort if the number of processes given by user do not match with the dat file
    if (procData.numRanks != procsPerAxis[0]*procsPerAxis[1]*procsPerAxis[2]) {
        if (procData.rank == 0) {
            printf("ERROR: The number of processes assigned in the dat file, %d, "
            "do not match the number of processes given as input, %d, while running the code\n",procsPerAxis[0]*procsPerAxis[1]*procsPerAxis[2],procData.numRanks);
        }
        finaliseMPI(&procData);
        return -1;
    }

    // INFO printing
    if(procData.rank == 0)
    {
#ifdef NDEBUG
    	printf("INFO: The compiler directive NDEBUG is enabled. Faster execution time is gained, "
    			"at the cost of less correctness checks during runtime!\n");
#else
		printf("INFO: The compiler directive NDEBUG is disabled. Checks for "
				"correctness are carried out at the cost of execution speed!\n"
				"      Use \"make speed\" for a faster execution time.\n");
#endif
    }

    // Domain decomposition & setting up neighbors
    domainDecompositionAndNeighbors(&procData, xlength, procsPerAxis);

    /*Allocate memory to pointers*/
    int totalsize = (procData.xLength[0]+2)*(procData.xLength[1]+2)*(procData.xLength[2]+2);
    c[0].collideField  = (double *)  malloc(Q*totalsize * sizeof( double ));
    c[0].streamField   = (double *)  malloc(Q*totalsize * sizeof( double ));

    /* calloc: only required to set boundary values. Sets every value to zero*/
    flagField     = (int *)  calloc(totalsize, sizeof( int ));

#ifdef UNITTEST
    double *massBefore[g_numComp];
    double *massAfter[g_numComp];
    double momentumBefore[3];
    double momentumAfter[3];

    for (int i = 0; i < g_numComp; i++) {
        massBefore[i] = (double *)  malloc(totalsize * sizeof( double ));
        massAfter[i] = (double *)  malloc(totalsize * sizeof( double ));
    }
#endif

    // Initialise the fields for all components
    initialiseComponents(c, g_numComp, flagField, &procData);

    // Allocate memory to send and read buffers
    initialiseBuffers(sendBuffer, readBuffer, procData.xLength, procData.neighbours, procData.bufferSize);

    // INFO Printing
    if(procData.rank == 0)
        printf("\nINFO: Storing cell data in VTS files.\n      Please use the"
        " \"Cell Data to Point Data\" filter in paraview to view nicely interpolated data. \n\n");

    // // Write the VTK at t = 0
    printf("R %i INFO: write vts file at time t = %d \n", procData.rank, t);
    writeVtsOutput(c, g_numComp, flagField, fName, t, xlength, &procData, procsPerAxis);

    // Combine VTS file at t = 0
    // Only done by root
    if (procData.rank == 0) {
        p_writeCombinedPVTSFile(g_numComp, fName, t, xlength, procsPerAxis);
    }

    beginProcTime = MPI_Wtime();

    for(t = 1; t <= timesteps; t++){

        //do extraction , swap , injection for - left/right, top/bottom, front/back for each component
        communicateComponents(sendBuffer, readBuffer, c, g_numComp, &procData);

        // Perform local streaming for each component
	    streamComponents(c, g_numComp, flagField, procData.xLength);

        // Swap the local fields for each component
        swapComponentFields(c, g_numComp);

#ifdef UNITTEST
        storeMassVector(c, g_numComp, massBefore, procData.xLength);
        computeGlobalMomentum(c, g_numComp, procData.xLength, momentumBefore);
#endif

        // Perform local collision
	    doCollision(c, procData.xLength);

#ifdef UNITTEST
        storeMassVector(c, g_numComp, massAfter, procData.xLength);
        computeGlobalMomentum(c, g_numComp, procData.xLength, momentumAfter);

        checkMassVector(massBefore, massAfter, procData.xLength, g_numComp, procData.rank);
        if (procData.rank == 0) {
            checkMomentum(momentumBefore, momentumAfter, g_numComp);
        }
#endif

        // Treat local boundaries for each component
	    treatComponentBoundary(c, g_numComp, flagField, &procData, sendBuffer, readBuffer);

        // Print VTS files at given interval
	    if (t%timestepsPerPlotting == 0){
            printf("R %i, INFO: write vts file at time t = %d \n", procData.rank, t);
	        writeVtsOutput(c, g_numComp, flagField, fName, t, xlength, &procData, procsPerAxis);
            // Combine VTS file at t
            if (procData.rank == 0) {
                p_writeCombinedPVTSFile(g_numComp, fName, t, xlength, procsPerAxis);
            }
	    }
    }
    endProcTime = MPI_Wtime();

    //get earliest start and latest finish of processors
    if(procData.rank == 0)
        printf("\nINFO: Please open the WS4_combined.*.pvts file to view the combined result.\n");

    MPI_Reduce(&beginProcTime, &beginSimTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&endProcTime, &endSimTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(procData.rank == 0){
    	double elapsedTime = endSimTime - beginSimTime;
    	int domTotalsize = (xlength+2)*(xlength+2)*(xlength+2);
    	printf("\n===============================================================\n");
		printf("\nINFO TIMING:\n");
		printf("Execution time (main loop): \t\t %.3f seconds \n", elapsedTime);
		printf("#cells (including boundary): \t\t %i cells \n", domTotalsize);
		printf("Mega Lattice Updates per Seconds: \t %f MLUPS \n",
				(domTotalsize/(1000000*elapsedTime))*timesteps);
		printf("\n===============================================================\n");
    }

    //free allocated heap memory
    free(c[0].streamField);
    free(c[0].collideField);
    free(flagField);

#ifdef UNITTEST
    for (int i = 0; i < g_numComp; i++) {
        free(massBefore[i]);
        free(massAfter[i]);
    }
#endif

    for (int i = LEFT; i <= BACK; i++) {
    	//set to NULL in initializeBuffers if buffer is not allocated (due to non existing neighbour)
        if(sendBuffer[i] != NULL) free(sendBuffer[i]);
        if(readBuffer[i] != NULL) free(readBuffer[i]);
    }

    finaliseMPI(&procData);
    return 0;
}

#endif
