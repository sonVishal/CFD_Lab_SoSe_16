#ifndef _MAIN_C_
#define _MAIN_C_
#include "streamCollide.h"
#include "collision.h"
#include "initLB.h"
#include "streaming.h"
#include "visualLB.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "common.h"
#include "parallel.h"
#include "unitTest.h"
#include <stdio.h>
#include <mpi/mpi.h>


//--------------------------------------------------------
//                        NOTE
//--------------------------------------------------------
// The framework is set up for doing multiple components, 
// but in the current state that is not implemented 
// properly.
//--------------------------------------------------------

//TODO: (TKS) We assume that mass is 1 and use numDensity as density in different places.
//            This must be fixed if we are doing several components.
//TODO: (TKS) Rename streamField and collideField.
//TODO: (TKS) Read in parameters for components from file.
//TODO: (TKS) Comment the code.
// TODO:(VS)   Remove velocity and reynolds number stuff after testing
// TODO: (TKS) Go through all ifndef/ifdef.

int main(int argc, char *argv[]){

    // Simulation parameters
    int xlength;
    int t = 0;
    int timesteps;
    int timestepsPerPlotting;

    t_procData procData;

    int procsPerAxis[3];

    // Send and read buffers for all possible directions:
    // Look at enum for index and direction correlation
    double *sendBuffer[6];
    double *readBuffer[6];

    initialiseMPI(&procData.rank,&procData.numRanks,argc,argv);

    // File printing parameters
    char fName[80];
    snprintf(fName, 80, "pv_files/project_rank");

    //Timing variables:
    double beginProcTime, endProcTime, beginSimTime, endSimTime;

    if (procData.rank == 0) {
        readNumComp(argc, &argv[1]);
    }else{
        MPI_Bcast(&numComp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    // Array of the components in our simulation
    t_component c[numComp];

    // Interacting potential
    double G[numComp][numComp];


    /* Read parameters*/
    //Only performed by the root and then broadcasted in 'broadcastValues'
    if (procData.rank == 0) {
        readParameters(&xlength, c, G, procsPerAxis, &timesteps, &timestepsPerPlotting, argc, &argv[1]);
    }

    // Broadcast the data from rank 0 (root) to other processes
    //TODO: (DL) maybe initialize also procData in here (collideField & streamField) & validity tests
    broadcastValues(procData.rank, &xlength, c, G, procsPerAxis, &timesteps, &timestepsPerPlotting);

    // Abort if the number of processes given by user do not match with the dat file
    //TODO: (DL) could we get this somewhere else than main?
    if (procData.numRanks != procsPerAxis[0]*procsPerAxis[1]*procsPerAxis[2]) {
        if (procData.rank == 0) {
            printf("ERROR: The number of processes assigned in the dat file, %d, "
            "do not match the number of processes given as input, %d, while running the code\n"
            ,procsPerAxis[0]*procsPerAxis[1]*procsPerAxis[2],procData.numRanks);
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

        // INFO Printing
        printf("\nINFO: Storing cell data in VTS files.\n      Please use the"
        " \"Cell Data to Point Data\" filter in paraview to view nicely interpolated data. \n\n");
    }

    // Domain decomposition & setting up neighbors
    domainDecompositionAndNeighbors(&procData, xlength, procsPerAxis);

    // TODO: Make a note about this not being a good way and we should calculate velocity etc?
    /*Allocate memory to pointers*/
    
    //------------------------------------------------------
    //          NOTE
    //------------------------------------------------------
    // We are aware that the current scheme is not optimal.
    // feq and the elements needed to calculate it should
    // be calculated on each process instead of it being
    // communicated and velocity etc. stored.
    //
    // The reason this is done is because of a very subtle
    // bug in our earlier code that we could not find, even
    // though we spent a week looking for it. We then chose
    // to restructure out code to be sure everything was
    // working as it should and unfortunately did not have
    // time for optimizations. Since the goal was to 
    // implement multiphase we feel this is fine enough
    //------------------------------------------------------

    int totalsize = (procData.xLength[0]+2)*(procData.xLength[1]+2)*(procData.xLength[2]+2);
    for (int i = 0; i < numComp; i++) {
        c[i].collideField = (double *)   malloc(Q*totalsize * sizeof( double ));
        c[i].streamField  = (double *)   malloc(Q*totalsize * sizeof( double ));
        c[i].feq          = (double *)   malloc(Q*totalsize * sizeof( double ));
        c[i].rho          = (double *)   malloc(totalsize   * sizeof( double ));
        c[i].velocity     = (double **)  malloc(totalsize   * sizeof( double* ));
        c[i].force        = (double **)  malloc(totalsize   * sizeof( double* ));

        for(int j = 0; j < totalsize; ++j){
            c[i].velocity[j] = (double*) malloc(3* sizeof( double));
            c[i].force[j]    = (double*) malloc(3* sizeof( double));
        }
    }

    //TODO: Make unit tests (?)
    //initializeUnitTest(totalsize);


    //----------------------------------------------------
    //                  NOTE
    //----------------------------------------------------
    // Currently flagField is not used in the simulation, 
    // but kept as a future extension to incorporate 
    // other boundary conditions
    //----------------------------------------------------

    /* calloc: only required to set boundary values. Sets every value to zero*/
    //int *flagField = (int *) calloc(totalsize, sizeof( int ));
    int *flagField = NULL;

    // Initialize all fields
    initialiseProblem(c, flagField, &procData);

    // #ifndef NDEBUG
    // initializeUnitTest(totalsize);
    // #endif

    // Allocate memory to send and read buffers
    //TODO: (DL) maybe hand in only procData for consistency?
    initialiseBuffers(sendBuffer, readBuffer, procData.xLength, procData.neighbours, procData.bufferSize);

    //Write the VTK at t = 0
    printf("R %i INFO: write vts file at time t = %d \n", procData.rank, t);
    writeVtsOutput(c, fName, t, xlength, &procData, procsPerAxis);

    // Combine VTS file at t = 0  -- only done by root
    if (procData.rank == 0) {
        p_writeCombinedPVTSFile(fName, t, xlength, procsPerAxis);
    }

    beginProcTime = MPI_Wtime();

    // Main time loop
    for(t = 1; t <= timesteps; t++){

        // Communicate required fields components to parallel boundaries 
        // and periodic boundaries
        communicateComponents(sendBuffer, readBuffer, c, &procData);

        // Collide and stream step
        streamCollide(c, &procData);

        // Compute the force between cells
		computeForce(c, &procData, flagField, G);

        // Compute density and velocity
        computeDensityAndVelocity(c, &procData);

        // Compute the equilibrium distributions
        computeFeq(c, &procData);

        // Swap fields. new --> old
        double *swap = NULL;
        for (int k = 0; k < numComp; ++k) {
            swap              = c[k].streamField;
            c[k].streamField  = c[k].collideField;
            c[k].collideField = swap;
        }

        
        // Print VTS files at given interval
		if (t%timestepsPerPlotting == 0){

            // TODO: Write something about this in input file, maybe
            //       make it a setting.
            //Double timestep for write in each iteration
            timestepsPerPlotting = timestepsPerPlotting*2;

            printf("R %i, INFO: write vts file at time t = %d \n", procData.rank, t);
            writeVtsOutput(c, fName, t, xlength, &procData, procsPerAxis);
            
            // Combine VTS file at t
            if (procData.rank == 0) {
                p_writeCombinedPVTSFile(fName, t, xlength, procsPerAxis);
            }
	    }

        if(t%10 == 0 && procData.rank == 0)
            printf("R %d Time t = %d done\n",procData.rank, t);
    }
    endProcTime = MPI_Wtime();

    if(procData.rank == 0)
        printf("\nINFO: Please open the WS4_combined.*.pvts file to view the combined result.\n");

    //Get earliest start and latest finish of processors
    MPI_Reduce(&beginProcTime, &beginSimTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&endProcTime,   &endSimTime,   1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(procData.rank == 0){
    	double elapsedTime = endSimTime - beginSimTime;
    	int domTotalsize   = (xlength+2)*(xlength+2)*(xlength+2);
    	printf("\n===============================================================\n");
		printf("\nINFO TIMING:\n");
		printf("Execution time (main loop): \t\t %.3f seconds \n", elapsedTime);
		printf("#cells (including boundary): \t\t %i cells \n", domTotalsize);
		printf("Mega Lattice Updates per Seconds: \t %f MLUPS \n",
				(domTotalsize/(1000000*elapsedTime))*timesteps);
		printf("\n===============================================================\n");
    }

    //free allocated heap memory
    //TODO: Remember to free memory
    for (int i = 0; i < numComp; i++) {
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

    for (int i = LEFT; i <= BACK; i++) {
    	//set to NULL in initializeBuffers if buffer is not allocated (due to non existing neighbour)
        if(sendBuffer[i] != NULL) free(sendBuffer[i]);
        if(readBuffer[i] != NULL) free(readBuffer[i]);
    }

    // #ifndef NDEBUG
    // freeUnitTest();
    // #endif

    finaliseMPI(&procData);
    return 0;
}

#endif
