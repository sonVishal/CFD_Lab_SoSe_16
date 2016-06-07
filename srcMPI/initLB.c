#include "helper.h"
#include "initLB.h"
#include "boundary.h"
#include <mpi/mpi.h>
#include <math.h>

int readParameters(int *xlength, double *tau, double *velocityWall, int *procsPerAxis,
	int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){

	if(argc != 2){
		char msg[200];
		snprintf(msg, 200, "There are %i arguments provided to the simulation. Only 1 argument "
				"providing \nthe path+filename (of relative to working directory) of the scenario "
				"is accepted!", argc-1);
		ERROR(msg);
	}

    double xvelocity, yvelocity, zvelocity;
	int iProc, jProc, kProc;
    double Re, u_wall, machNr;

    /* Read values from file given in argv */
    READ_INT(*argv, *xlength);
    READ_DOUBLE(*argv, Re);

    READ_DOUBLE(*argv, xvelocity);
    READ_DOUBLE(*argv, yvelocity);
    READ_DOUBLE(*argv, zvelocity);

    velocityWall[0] = xvelocity;
    velocityWall[1] = yvelocity;
    velocityWall[2] = zvelocity;

    READ_INT(*argv, *timesteps);
    READ_INT(*argv, *timestepsPerPlotting);

	READ_INT(*argv, iProc);
	READ_INT(*argv, jProc);
	READ_INT(*argv, kProc);

	procsPerAxis[0] = iProc;
	procsPerAxis[1] = jProc;
	procsPerAxis[2] = kProc;

    /*Calculates tau from the Reynolds number*/
    u_wall  = sqrt(xvelocity*xvelocity + yvelocity*yvelocity+zvelocity*zvelocity);
    *tau    =  u_wall*(*xlength)/(C_S*C_S*Re)+0.5;
    machNr  = u_wall/C_S;

    printf("\nINFO: Calculated tau = %f \n", *tau);
    printf("\nINFO: Wall speed = %f \n", u_wall);
    printf("\nINFO: Mach number = %f \n\n", machNr);

    /* valid settings check*/
    if(u_wall >= C_S){
    	ERROR("Wall speed is supersonic (aborting)! \n");
    }

    if(*tau<=0.5 || *tau>2){
        ERROR("Tau is out of stability region (0.5,2.0) (aborting)! \n");
    }

    /*We allow user defined mach number tolerance for Ma << 1 (default = 0.1)
      To change please look at LBDefinitions.h*/
    if(machNr > machNrTol){
        char buffer[80];
        snprintf(buffer, 80, "Mach number is larger than %f (aborting)! \n",machNrTol);
    	ERROR(buffer);
    }

	// Check if iProc, jProc and kProc are less than the domain size and nonzero
	if (iProc <= 0 || iProc > *xlength || jProc <= 0 || jProc > *xlength ||
		kProc <= 0 || kProc > *xlength) {
		char buffer[80];
        snprintf(buffer, 80, "Please make sure that iProc, jProc and kProc are greater than 0 and less than xlength = %d (aborting)! \n",*xlength);
    	ERROR(buffer);
	}

	if (*xlength%iProc != 0 || *xlength%jProc != 0 || *xlength%kProc != 0) {
		char buffer[160];
        snprintf(buffer, 160, "WARNING: The domain decomposition is not uniform. This might create load unbalances between processes.\n");
	}

  return 0;
}

void initialiseFields(double *collideField, double *streamField, int *flagField, t_procData thisProcData){

    /*Setting initial distributions*/
    //f_i(x,0) = f^eq(1,0,0) = w_i

    // current cell index
    int idx;

    // Temporary variables for xlength^2
    int const xylen = (thisProcData.xLength[0]+2)*(thisProcData.xLength[1]+2);

    // Temporary variables for z and y offsets
    int zOffset, yzOffset;

    /* initialize collideField and streamField */
    int x,y,z;
    for ( z = 0; z <= thisProcData.xLength[2]+1; ++z) {
        zOffset = z*xylen;
        for ( y = 0; y <= thisProcData.xLength[1]+1; ++y) {
            yzOffset = y*(thisProcData.xLength[0]+2) + zOffset;
            for ( x = 0; x <= thisProcData.xLength[0]+1; ++x) {
                // Compute the base index
                idx = Q*(yzOffset + x);
                for (int i = 0; i < Q; ++i) {
                    collideField[idx+i] = LATTICEWEIGHTS[i];
                    streamField[idx+i]  = LATTICEWEIGHTS[i];
                }
            }
        }
    }


    /*Looping over boundary of flagFields*/
    //All points set to zero at memory allocation (using calloc)
	int endOuter, endInner, fixedValue, boundaryType, wallIdx;
	// First set up the ghost boundary layer
	for (wallIdx = LEFT; wallIdx <= BACK; wallIdx++) {
		if (thisProcData.neighbours[wallIdx] == MPI_PROC_NULL) {
			p_setIterationParameters(&endOuter, &endInner, &fixedValue, thisProcData, wallIdx);
			if (wallIdx == TOP) {
				boundaryType = MOVING_WALL;
			} else {
				boundaryType = NO_SLIP;
			}
			for (z = 0; z <= endOuter; z++) {
				for (y = 0; y <= endInner; y++) {
					idx = p_computeCellOffset(z,y,fixedValue,thisProcData.xLength,wallIdx);
					flagField[idx] = boundaryType;
				}
			}
		}
	}

	// Now set up the parallel boundary layer
	for (wallIdx = LEFT; wallIdx <= BACK; wallIdx++) {
		if (thisProcData.neighbours[wallIdx] != MPI_PROC_NULL) {
			p_setIterationParameters(&endOuter, &endInner, &fixedValue, thisProcData, wallIdx);
			for (z = 0; z <= endOuter; z++) {
				for (y = 0; y <= endInner; y++) {
					idx = p_computeCellOffset(z,y,fixedValue,thisProcData.xLength,wallIdx);
					flagField[idx] = PARALLEL_BOUNDARY;
				}
			}
		}
	}
}

void initialiseMPI(int *rank, int *numRanks, int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,numRanks);
    MPI_Comm_rank(MPI_COMM_WORLD,rank);
}

void p_domainDecompositionAndNeighbors(t_procData *procData, const int xlength, const int * const procsPerAxis) {
	// printf("Calling domain decomposition from %d\n",procData->rank);
	int procPos[3] = {0,0,0};
	p_rankToPos(procsPerAxis,procData->rank,procPos);

    /* Compute the subdomain size and save it into procData */
 	procData->xLength[0] = xlength/procsPerAxis[0];
    procData->xLength[1] = xlength/procsPerAxis[1];
    procData->xLength[2] = xlength/procsPerAxis[2];
    // If the proc is at the end of some axis then add the remaining length
    procData->xLength[0] += (procPos[0] == procsPerAxis[0]-1)?xlength%procsPerAxis[0]:0;
    procData->xLength[1] += (procPos[1] == procsPerAxis[1]-1)?xlength%procsPerAxis[1]:0;
    procData->xLength[2] += (procPos[2] == procsPerAxis[2]-1)?xlength%procsPerAxis[2]:0;


    /* Decide whether it is a ghost boundary (= MPI_PROC_NULL) or a parallel boundary (rank of neighbour) */
    procData->neighbours[LEFT]   = (procPos[1] == 0) 				  ? MPI_PROC_NULL : procData->rank-procsPerAxis[0];
    procData->neighbours[RIGHT]  = (procPos[1] == procsPerAxis[1]-1)  ? MPI_PROC_NULL : procData->rank+procsPerAxis[0];

    procData->neighbours[BACK]   = (procPos[0] == 0) 				  ? MPI_PROC_NULL : procData->rank-1;
    procData->neighbours[FRONT]  = (procPos[0] == procsPerAxis[0]-1)  ? MPI_PROC_NULL : procData->rank+1;

    procData->neighbours[BOTTOM] = (procPos[2] == 0) 				  ? MPI_PROC_NULL : procData->rank-procsPerAxis[1]*procsPerAxis[0];
    procData->neighbours[TOP]    = (procPos[2] == procsPerAxis[2]-1)  ? MPI_PROC_NULL : procData->rank+procsPerAxis[1]*procsPerAxis[0];

//    if (procPos[0] == 0) {
//        procData->neighbours[BACK] = MPI_PROC_NULL;
//    } else {
//        procData->neighbours[BACK] = procData->rank-1;
//    }
//    // Front
//    if (procPos[0] == procsPerAxis[0]-1) {
//        procData->neighbours[FRONT] = MPI_PROC_NULL;
//    } else {
//        procData->neighbours[FRONT] = procData->rank+1;
//    }
//    // Left
//    if (procPos[1] == 0) {
//        procData->neighbours[LEFT] = MPI_PROC_NULL;
//    } else {
//        procData->neighbours[LEFT] = procData->rank-procsPerAxis[0];
//    }
//    // Right
//    if (procPos[1] == procsPerAxis[1]-1) {
//        procData->neighbours[RIGHT] = MPI_PROC_NULL;
//    } else {
//        procData->neighbours[RIGHT] = procData->rank+procsPerAxis[0];
//    }
//    // Bottom
//    if (procPos[2] == 0) {
//        procData->neighbours[BOTTOM] = MPI_PROC_NULL;
//    } else {
//        procData->neighbours[BOTTOM] = procData->rank-procsPerAxis[1]*procsPerAxis[0];
//    }
//    // Top
//    if (procPos[2] == procsPerAxis[2]-1) {
//        procData->neighbours[TOP] = MPI_PROC_NULL;
//    } else {
//        procData->neighbours[TOP] = procData->rank+procsPerAxis[1]*procsPerAxis[0];
//    }
	// printf("Proc %d \t Position (%d,%d,%d)\n",procData->rank,procPos[0],procPos[1],procPos[2]);
	// printProcDataPos(*procData,procPos);
}

void initialiseBuffers(double *sendBuffer[6], double *readBuffer[6], int *xlength) {

    //TODO: (TKS) Consider allocating only the amount of elements needed.
	int bufferSize		= 5*(xlength[0]*xlength[2])*sizeof(double);
	sendBuffer[LEFT] 	= (double *) malloc(bufferSize);
	sendBuffer[RIGHT] 	= (double *) malloc(bufferSize);
	readBuffer[LEFT] 	= (double *) malloc(bufferSize);
	readBuffer[RIGHT] 	= (double *) malloc(bufferSize);


	bufferSize 			= 5*((xlength[0]+2)*xlength[1])*sizeof(double);
	sendBuffer[TOP] 	= (double *) malloc(bufferSize);
	sendBuffer[BOTTOM] 	= (double *) malloc(bufferSize);
	readBuffer[TOP] 	= (double *) malloc(bufferSize);
	readBuffer[BOTTOM] 	= (double *) malloc(bufferSize);

	bufferSize 			= 5*((xlength[1]+2)*(xlength[2]+2))*sizeof(double);
	sendBuffer[FRONT] 	= (double *) malloc(bufferSize);
	sendBuffer[BACK] 	= (double *) malloc(bufferSize);
	readBuffer[FRONT] 	= (double *) malloc(bufferSize);
	readBuffer[BACK] 	= (double *) malloc(bufferSize);
}

/* TODO: (DL) This should go somewhere else!! */
void finaliseMPI(t_procData *procData) {
    fflush(stdout);
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Rank %i: FINISHED. \n", procData->rank);
    MPI_Finalize();
}
