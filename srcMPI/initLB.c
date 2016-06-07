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

	for (x = LEFT; x <= BACK; x++) {
		if (thisProcData.neighbours[x] == MPI_PROC_NULL) {

		}
	}

    //These are the no-slip walls
    //fixed: z = 0

    for (y = 0; y <= thisProcData.xLength[1]+1; y++) {
        idx = y*(thisProcData.xLength[0]+2);
        for (x = 0; x <= thisProcData.xLength[1]+1; x++) {
            flagField[x+idx] = 1;
        }
    }

    //fixed: x = 0
    //We start at 1 to not include previous cells again from z = 0
    for (z = 1; z <= thisProcData.xLength[2]; z++) {
        zOffset = z*xylen;
        for (y = 0; y <= thisProcData.xLength[1]+1; y++) {
            flagField[zOffset+y*(thisProcData.xLength[0]+2)] = 1;
        }
    }

    //fixed: x = thisProcData.xLength+1
    //We start at 1 to not include previous cells again from z = 0
    for (z = 1; z <= thisProcData.xLength[2]; z++) {
        zOffset = z*xylen + thisProcData.xLength[0] + 1;
        for (y = 0; y <= thisProcData.xLength[1]+1; y++) {
            flagField[zOffset+y*(thisProcData.xLength[0]+2)] = 1;
        }
    }

    //fixed: y = 0
    //from 1:thisProcData.xLength only, to not include cells at upper, lower, left and right edges
    //The edge cells are set in the other loops.
    for (z = 1; z <= thisProcData.xLength[2]; z++) {
        zOffset = z*xylen;
        for (x = 1; x <= thisProcData.xLength[0]; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    //fixed: y = thisProcData.xLength+1
    //same reasoning for index range as in fixed y=0
    for (z = 1; z <= thisProcData.xLength[2]; z++) {
        zOffset = z*xylen + (thisProcData.xLength[1]+1)*(thisProcData.xLength[0]+2);
        for (x = 1; x <= thisProcData.xLength[0]; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    // This is the moving wall. All cells at z=thisProcData.xLength+1 are included (also the edge cells).
    // fixed: z = thisProcData.xLength+1
    zOffset = (thisProcData.xLength[2]+1)*xylen;
    for (y = 0; y <= thisProcData.xLength[1]+1; y++) {
        idx = zOffset + y*(thisProcData.xLength[0]+2);
        for (x = 0; x <= thisProcData.xLength[0]+1; x++) {
            flagField[x+idx] = 2;
        }
    }
}

void initialiseMPI(int *rank, int *numRanks, int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,numRanks);
    MPI_Comm_rank(MPI_COMM_WORLD,rank);
}

void p_domainDecompositionAndNeighbors(t_procData *procData, const int xlength, const int * const procsPerAxis) {
	int procPos[3] = {0,0,0};
	procData->xLength[0] = xlength/procsPerAxis[0];
    procData->xLength[1] = xlength/procsPerAxis[1];
    procData->xLength[2] = xlength/procsPerAxis[2];
    p_indexToPos(procsPerAxis,procData->rank,procPos);
    // If the proc is at the end of some axis then add the remaining length
    procData->xLength[0] += (procPos[0] == procsPerAxis[0]-1)?xlength%procsPerAxis[0]:0;
    procData->xLength[1] += (procPos[1] == procsPerAxis[1]-1)?xlength%procsPerAxis[1]:0;
    procData->xLength[2] += (procPos[2] == procsPerAxis[2]-1)?xlength%procsPerAxis[2]:0;
	// Back
    if (procPos[0] == 0) {
        procData->neighbours[BACK] = MPI_PROC_NULL;
    } else {
        procData->neighbours[BACK] = procData->rank-1;
    }
    // Front
    if (procPos[0] == procsPerAxis[0]-1) {
        procData->neighbours[FRONT] = MPI_PROC_NULL;
    } else {
        procData->neighbours[FRONT] = procData->rank+1;
    }
    // Left
    if (procPos[1] == 0) {
        procData->neighbours[LEFT] = MPI_PROC_NULL;
    } else {
        procData->neighbours[LEFT] = procData->rank-procsPerAxis[0];
    }
    // Right
    if (procPos[1] == procsPerAxis[1]-1) {
        procData->neighbours[RIGHT] = MPI_PROC_NULL;
    } else {
        procData->neighbours[RIGHT] = procData->rank+procsPerAxis[0];
    }
    // Bottom
    if (procPos[2] == 0) {
        procData->neighbours[BOTTOM] = MPI_PROC_NULL;
    } else {
        procData->neighbours[BOTTOM] = procData->rank+procsPerAxis[1]*procsPerAxis[0];
    }
    // Top
    if (procPos[2] == procsPerAxis[2]-1) {
        procData->neighbours[TOP] = MPI_PROC_NULL;
    } else {
        procData->neighbours[TOP] = procData->rank-procsPerAxis[1]*procsPerAxis[0];
    }
}

void initialiseBuffers(double *sendBuffer[6], double *readBuffer[6], int *xlength) {
	ERROR("TODO");
}

void finaliseMPI() {
	ERROR("TODO");
}
