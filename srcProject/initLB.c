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

int readParameters(int *xlength, t_component *c, double G[numComp][numComp], int *procsPerAxis,
    int *timesteps, int *timestepsPerPlotting,
	int argc, char *argv[]){

	if(argc != 2){
		char msg[200];
		snprintf(msg, 200, "There are %i arguments provided to the simulation. Only 1 argument "
				"providing \nthe path+filename (of relative to working directory) of the scenario "
				"is accepted!", argc-1);
		ERROR(msg);
	}

	int iProc, jProc, kProc;

    // [> Read values from file given in argv <]
    READ_INT(*argv, *xlength);

	// Read the tau and mass for each component
	// * If mass and tau are provided for less number of components
	//   read_double will give an error
	// * If mass and tau are provided for more number of components
	// 	 they will not be read due to the for loop
	for (int i = 0; i < numComp; i++) {
		char tempName[20];
		snprintf(tempName, 20, "tau%d",i);
		read_double(*argv, tempName, &c[i].tau);

		snprintf(tempName, 20, "m%d",i);
		read_double(*argv, tempName, &c[i].m);

		for (int j = 0; j < numComp; j++) {
			snprintf(tempName, 20, "G%d%d",i,j);
			read_double(*argv, tempName, &G[i][j]);
		}

		snprintf(tempName, 20, "psi%d", i);
		read_int(*argv, tempName, &c[i].psiFctCode);
	}

    READ_INT(*argv, *timesteps);
    READ_INT(*argv, *timestepsPerPlotting);

	READ_INT(*argv, iProc);
	READ_INT(*argv, jProc);
	READ_INT(*argv, kProc);

	procsPerAxis[0] = iProc;
	procsPerAxis[1] = jProc;
	procsPerAxis[2] = kProc;

    //TODO: make tau check again (can be multiple...)
    // if(*tau<=0.5 || *tau>2){
    //     ERROR("Tau is out of stability region (0.5,2.0) (aborting)! \n");
    // }

    // [>We allow user defined mach number tolerance for Ma << 1 (default = 0.1)
    //   To change please look at LBDefinitions.h*/
    // if(machNr > machNrTol){
    //     char buffer[80];
    //     snprintf(buffer, 80, "Mach number is larger than %f (aborting)! \n",machNrTol);
    //     ERROR(buffer);
    // }

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


// NOTE: Currently only supports random initialization for multiphase
void initialiseFields(t_component * c, const t_procData * const procData){

    // [>Setting initial distributions<]
    //f_i(x,0) = f^eq(random,0)

    // current cell index
    int fieldIdx, cellIdx;

	// How much initial difference is allowed in density
	double rhoVar = 0.01*rhoRef;  //Shan Chen

	// Initially velocity is 0
	double u0[3] = {0.0, 0.0, 0.0};

	//TODO: (DL) this 'if' is for reference solution, can be deleted when cleaning up...
	if(procData->rank == 0){
		srand(1);
	}else if(procData->rank == 1){
		srand(3);
	}else{
		srand(procData->rank*10);
	}

    // [> initialize collideField and streamField <]
    int x,y,z;
    for ( z = 1; z <= procData->xLength[2]; ++z) {
        for ( y = 1; y <= procData->xLength[1]; ++y) {
            for ( x = 1; x <= procData->xLength[0]; ++x) {
                // Compute the base index
				fieldIdx = p_computeCellOffsetXYZ(x, y, z, procData->xLength);
                cellIdx = Q*fieldIdx;

				double rnd = ((double)rand()/(double)RAND_MAX);
				c->rho[fieldIdx] = rhoRef - 0.5*rhoVar + rhoVar*rnd; //Shan Chen
				//c->rho[fieldIdx] = rhoRef + rnd; //Sukop

				if(x == 1 && y == 1 && z == 2 && procData->rank == 1){
					printf("R%i: @(1,1,2) %.16f \n", procData->rank, c->rho[fieldIdx]);
				}

				if(x == 1 && y == 1 && z == 1 && procData->rank == 0){
					printf("R%i: @(1,1,1) %.16f \n", procData->rank, c->rho[fieldIdx]);
				}


				computeFeqCell(&(c->rho[fieldIdx]), u0, &(c->feq[cellIdx]));
                for (int i = 0; i < Q; ++i) {
                    c->collideField[cellIdx+i] = c->feq[cellIdx+i];
                    c->streamField[cellIdx+i]  = c->feq[cellIdx+i];
                }
            }
        }
    }
}

// Currently the initialization is only for single component multiphase simulation
// TODO: Box and random initialization for multicomponents as per Sukop code
void initialiseComponents(t_component *c, int *flagField, const t_procData * const procData) {

	// Initialize the collide and stream fields for each component
	// TODO: this also changes in multicomponent
	// Set up a seed based on the procData
	// srand(procData->rank*(procData->rank+1)+1);
	for (int i = 0; i < numComp; i++) {
		initialiseFields(&c[i], procData);
	}

	// Flag field is component independent

	// [>Looping over boundary of flagFields<]
    //All points set to zero at memory allocation (using calloc)
	int innerIt,outerIt;	// Iteration variables
	int endOuter, endInner, fixedValue, wallIdx, idx;

	// First set up the parallel boundary layer
	for (wallIdx = LEFT; wallIdx <= BACK; wallIdx++) {
		if (procData->neighbours[wallIdx] != MPI_PROC_NULL) {
			p_setIterationParameters(&endOuter, &endInner, &fixedValue, procData, wallIdx);
			for (outerIt = 0; outerIt <= endOuter; outerIt++) {
				for (innerIt = 0; innerIt <= endInner; innerIt++) {
					idx = p_computeCellOffset(outerIt,innerIt,fixedValue,procData->xLength,wallIdx);
					flagField[idx] = PARALLEL_BOUNDARY;
				}
			}
		}
	}

	// Now set up the ghost boundary layer with no slip
	for (wallIdx = LEFT; wallIdx <= BACK; wallIdx++) {
		if (procData->neighbours[wallIdx] == MPI_PROC_NULL && wallIdx != TOP) {
			p_setIterationParameters(&endOuter, &endInner, &fixedValue, procData, wallIdx);
			for (outerIt = 0; outerIt <= endOuter; outerIt++) {
				for (innerIt = 0; innerIt <= endInner; innerIt++) {
					idx = p_computeCellOffset(outerIt,innerIt,fixedValue,procData->xLength,wallIdx);
					flagField[idx] = PERIODIC_BOUNDARY;
				}
			}
		}
	}

	// Now set up the ghost boundary layer with moving wall - if it's at the TOP ghost layer
	if (procData->neighbours[TOP] == MPI_PROC_NULL) {
		p_setIterationParameters(&endOuter, &endInner, &fixedValue, procData, TOP);
		for (outerIt = 0; outerIt <= endOuter; outerIt++) {
			for (innerIt = 0; innerIt <= endInner; innerIt++) {
				idx = p_computeCellOffset(outerIt,innerIt,fixedValue,procData->xLength, TOP);
				flagField[idx] = PERIODIC_BOUNDARY;
			}
		}
	}
}
