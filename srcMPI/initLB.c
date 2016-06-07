#include "helper.h"
#include "initLB.h"
#include "LBDefinitions.h"
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

void initialiseFields(double *collideField, double *streamField, int *flagField,
	int *xlength, int rank, int numRanks){

    /*Setting initial distributions*/
    //f_i(x,0) = f^eq(1,0,0) = w_i

    // current cell index
    int idx;

    // Temporary variables for xlength^2
    int const xylen = (xlength[0]+2)*(xlength[1]+2);

    // Temporary variables for z and y offsets
    int zOffset, yzOffset;

    /* initialize collideField and streamField */
    int x,y,z;
    for ( z = 0; z <= xlength[2]+1; ++z) {
        zOffset = z*xylen;
        for ( y = 0; y <= xlength[1]+1; ++y) {
            yzOffset = y*(xlength[0]+2) + zOffset;
            for ( x = 0; x <= xlength[0]+1; ++x) {
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

    //These are the no-slip walls
    //fixed: z = 0
    for (y = 0; y <= xlength[1]+1; y++) {
        idx = y*(xlength[0]+2);
        for (x = 0; x <= xlength[1]+1; x++) {
            flagField[x+idx] = 1;
        }
    }

    //fixed: x = 0
    //We start at 1 to not include previous cells again from z = 0
    for (z = 1; z <= xlength[2]; z++) {
        zOffset = z*xylen;
        for (y = 0; y <= xlength[1]+1; y++) {
            flagField[zOffset+y*(xlength[0]+2)] = 1;
        }
    }

    //fixed: x = xlength+1
    //We start at 1 to not include previous cells again from z = 0
    for (z = 1; z <= xlength[2]; z++) {
        zOffset = z*xylen + xlength[0] + 1;
        for (y = 0; y <= xlength[1]+1; y++) {
            flagField[zOffset+y*(xlength[0]+2)] = 1;
        }
    }

    //fixed: y = 0
    //from 1:xlength only, to not include cells at upper, lower, left and right edges
    //The edge cells are set in the other loops.
    for (z = 1; z <= xlength[2]; z++) {
        zOffset = z*xylen;
        for (x = 1; x <= xlength[0]; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    //fixed: y = xlength+1
    //same reasoning for index range as in fixed y=0
    for (z = 1; z <= xlength[2]; z++) {
        zOffset = z*xylen + (xlength[1]+1)*(xlength[0]+2);
        for (x = 1; x <= xlength[0]; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    // This is the moving wall. All cells at z=xlength+1 are included (also the edge cells).
    // fixed: z = xlength+1
    zOffset = (xlength[2]+1)*xylen;
    for (y = 0; y <= xlength[1]+1; y++) {
        idx = zOffset + y*(xlength[0]+2);
        for (x = 0; x <= xlength[0]+1; x++) {
            flagField[x+idx] = 2;
        }
    }
}

void initialiseMPI(int *rank, int *numRanks, int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,numRanks);
    MPI_Comm_rank(MPI_COMM_WORLD,rank);
    ERROR("TODO");
}

void initialiseBuffers(double *sendBuffer[6], double *readBuffer[6], int *xlength) {
	ERROR("TODO");
}

void finaliseMPI() {
	ERROR("TODO");
}
