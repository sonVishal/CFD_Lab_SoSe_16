#include "initLB.h"

void readNumComp(int argc, char *argv[]) {
	if(argc != 2){
		char msg[200];
		snprintf(msg, 200, "There are %i arguments provided to the simulation. Only 1 argument "
				"providing \nthe path+filename (of relative to working directory) of the scenario "
				"is accepted!", argc-1);
		ERROR(msg);
	}

	READ_INT(*argv, numComp);

	MPI_Bcast(&numComp, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

// TODO: (VS) Remove velocity and reynolds number in final version from entire code
int readParameters(int *xlength, t_component *c, double G[numComp][numComp],
	double *velocityWall, int *procsPerAxis, int *timesteps, int *timestepsPerPlotting,
	int argc, char *argv[]){

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

	//Broadcast it immediately so that all the others have it
	MPI_Bcast(numComp, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Read the tau and mass for each component
	// * If mass and tau are provided for less number of components
	//   p_read_double will give an error
	// * If mass and tau are provided for more number of components
	// 	 they will not be read due to the for loop
	for (int i = 0; i < numComp; i++) {
		char tempName[20];
		snprintf(tempName, 20, "tau%d",i);
		p_read_double(*argv, tempName, &c[i].tau);
		snprintf(tempName, 20, "m%d",i);
		p_read_double(*argv, tempName, &c[i].m);
		for (int j = 0; j < numComp; j++) {
			snprintf(tempName, 20, "G%d%d",i,j);
			p_read_double(*argv, tempName, &G[i][j]);
		}
	}

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
    // *tau    =  u_wall*(*xlength)/(C_S*C_S*Re)+0.5;
    machNr  = u_wall/C_S;

    // printf("\nINFO: Calculated tau = %f \n", *tau);
    printf("\nINFO: Wall speed = %f \n", u_wall);
    printf("\nINFO: Mach number = %f \n\n", machNr);

    /* valid settings check*/
    if(u_wall >= C_S){
    	ERROR("Wall speed is supersonic (aborting)! \n");
    }

    // if(*tau<=0.5 || *tau>2){
    //     ERROR("Tau is out of stability region (0.5,2.0) (aborting)! \n");
    // }

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
		char buffer[120];
        snprintf(buffer, 120, "Please make sure that iProc, jProc and kProc are greater than 0 and less than xlength = %d (aborting)! \n",*xlength);
    	ERROR(buffer);
	}

	if (*xlength%iProc != 0 || *xlength%jProc != 0 || *xlength%kProc != 0) {
        printf("WARNING: The domain decomposition is not uniform. This might create load unbalances between processes.\n");
	}


  return 0;
}

void initialiseFields(double *collideField, double *streamField, const t_procData * const thisProcData){

    /*Setting initial distributions*/
    //f_i(x,0) = f^eq(1,0,0) = w_i

    // current cell index
    int idx;

    // Temporary variables for xlength^2
    int const xylen = (thisProcData->xLength[0]+2)*(thisProcData->xLength[1]+2);

    // Temporary variables for z and y offsets
    int zOffset, yzOffset;

    /* initialize collideField and streamField */
    int x,y,z;
    for ( z = 0; z <= thisProcData->xLength[2]+1; ++z) {
        zOffset = z*xylen;
        for ( y = 0; y <= thisProcData->xLength[1]+1; ++y) {
            yzOffset = y*(thisProcData->xLength[0]+2) + zOffset;
            for ( x = 0; x <= thisProcData->xLength[0]+1; ++x) {
                // Compute the base index
                idx = Q*(yzOffset + x);
                for (int i = 0; i < Q; ++i) {
                    collideField[idx+i] = LATTICEWEIGHTS[i];
                    streamField[idx+i]  = LATTICEWEIGHTS[i];

					// TODO: (VS) Remove.
					// This is good for debugging since we see different values
					// in each subdomain.
					// collideField[idx+i] = (double)(thisProcData->rank+1)/100.0;
                    // streamField[idx+i]  = (double)(thisProcData->rank+1)/100.0;
                }
            }
        }
    }
}

void initialiseComponents(t_component *c, int numComp, int *flagField, const t_procData * const thisProcData) {

	// Initialize the collide and stream fields for each component
	for (int i = 0; i < numComp; i++) {
		initialiseFields(c[i].collideField, c[i].streamField, thisProcData);
	}

	// Flag field is component independent

	/*Looping over boundary of flagFields*/
    //All points set to zero at memory allocation (using calloc)
	int innerIt,outerIt;	// Iteration variables
	int endOuter, endInner, fixedValue, wallIdx, idx;

	// First set up the parallel boundary layer
	for (wallIdx = LEFT; wallIdx <= BACK; wallIdx++) {
		if (thisProcData->neighbours[wallIdx] != MPI_PROC_NULL) {
			p_setIterationParameters(&endOuter, &endInner, &fixedValue, thisProcData, wallIdx);
			for (outerIt = 0; outerIt <= endOuter; outerIt++) {
				for (innerIt = 0; innerIt <= endInner; innerIt++) {
					idx = p_computeCellOffset(outerIt,innerIt,fixedValue,thisProcData->xLength,wallIdx);
					flagField[idx] = PARALLEL_BOUNDARY;
				}
			}
		}
	}

	//TODO: (DL) only quickly changed to PERIODIC_BC

	// Now set up the ghost boundary layer with no slip
	for (wallIdx = LEFT; wallIdx <= BACK; wallIdx++) {
		if (thisProcData->neighbours[wallIdx] == MPI_PROC_NULL && wallIdx != TOP) {
			p_setIterationParameters(&endOuter, &endInner, &fixedValue, thisProcData, wallIdx);
			for (outerIt = 0; outerIt <= endOuter; outerIt++) {
				for (innerIt = 0; innerIt <= endInner; innerIt++) {
					idx = p_computeCellOffset(outerIt,innerIt,fixedValue,thisProcData->xLength,wallIdx);
					flagField[idx] = PERIODIC_BOUNDARY;
				}
			}
		}
	}

	// Now set up the ghost boundary layer with moving wall - if it's at the TOP ghost layer
	if (thisProcData->neighbours[TOP] == MPI_PROC_NULL) {
		p_setIterationParameters(&endOuter, &endInner, &fixedValue, thisProcData, TOP);
		for (outerIt = 0; outerIt <= endOuter; outerIt++) {
			for (innerIt = 0; innerIt <= endInner; innerIt++) {
				idx = p_computeCellOffset(outerIt,innerIt,fixedValue,thisProcData->xLength, TOP);
				flagField[idx] = PERIODIC_BOUNDARY;
			}
		}
	}
}
