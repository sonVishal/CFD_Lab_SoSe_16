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
    double Re;

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

    READ_DOUBLE(*argv, Re);

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

//void initialiseFields(double *collideField, double *streamField, const t_procData * const thisProcData){

    //[>Setting initial distributions<]
    ////f_i(x,0) = f^eq(1,0,0) = w_i

    //// current cell index
    //int idx;

    //// Temporary variables for xlength^2
    //int const xylen = (thisProcData->xLength[0]+2)*(thisProcData->xLength[1]+2);

    //// Temporary variables for z and y offsets
    //int zOffset, yzOffset;

    //[> initialize collideField and streamField <]
    //int x,y,z;
    //for ( z = 0; z <= thisProcData->xLength[2]+1; ++z) {
        //zOffset = z*xylen;
        //for ( y = 0; y <= thisProcData->xLength[1]+1; ++y) {
            //yzOffset = y*(thisProcData->xLength[0]+2) + zOffset;
            //for ( x = 0; x <= thisProcData->xLength[0]+1; ++x) {
                //// Compute the base index
                //idx = Q*(yzOffset + x);
                //for (int i = 0; i < Q; ++i) {
                    //collideField[idx+i] = 200*LATTICEWEIGHTS[i] + (double)rand()/(double)RAND_MAX/Q;/[>(1+(double)rand()/(double)RAND_MAX/100);
                    //streamField[idx+i]  = 200*LATTICEWEIGHTS[i] + (double)rand()/(double)RAND_MAX/Q;/[>(1+(double)rand()/(double)RAND_MAX/100);

					//// TODO: (VS) Remove.
					//// This is good for debugging since we see different values
					//// in each subdomain.
					//// collideField[idx+i] = (double)(thisProcData->rank+1)/100.0;
                    //// streamField[idx+i]  = (double)(thisProcData->rank+1)/100.0;
                //}
            //}
        //}
    //}
//}

//void initialiseComponents(t_component *c, int *flagField, const t_procData * const thisProcData) {

	//// Initialize the collide and stream fields for each component
	//for (int i = 0; i < numComp; i++) {
		//initialiseFields(c[i].collideField, c[i].streamField, thisProcData);
	//}

	//// Flag field is component independent

	//[>Looping over boundary of flagFields<]
    ////All points set to zero at memory allocation (using calloc)
	//int innerIt,outerIt;	// Iteration variables
	//int endOuter, endInner, fixedValue, wallIdx, idx;

	//// First set up the parallel boundary layer
	//for (wallIdx = LEFT; wallIdx <= BACK; wallIdx++) {
		//if (thisProcData->neighbours[wallIdx] != MPI_PROC_NULL) {
			//p_setIterationParameters(&endOuter, &endInner, &fixedValue, thisProcData, wallIdx);
			//for (outerIt = 0; outerIt <= endOuter; outerIt++) {
				//for (innerIt = 0; innerIt <= endInner; innerIt++) {
					//idx = p_computeCellOffset(outerIt,innerIt,fixedValue,thisProcData->xLength,wallIdx);
					//flagField[idx] = PARALLEL_BOUNDARY;
				//}
			//}
		//}
	//}

	////TODO: (DL) only quickly changed to PERIODIC_BC

	//// Now set up the ghost boundary layer with no slip
	//for (wallIdx = LEFT; wallIdx <= BACK; wallIdx++) {
		//if (thisProcData->neighbours[wallIdx] == MPI_PROC_NULL && wallIdx != TOP) {
			//p_setIterationParameters(&endOuter, &endInner, &fixedValue, thisProcData, wallIdx);
			//for (outerIt = 0; outerIt <= endOuter; outerIt++) {
				//for (innerIt = 0; innerIt <= endInner; innerIt++) {
					//idx = p_computeCellOffset(outerIt,innerIt,fixedValue,thisProcData->xLength,wallIdx);
					//flagField[idx] = PERIODIC_BOUNDARY;
				//}
			//}
		//}
	//}

	//// Now set up the ghost boundary layer with moving wall - if it's at the TOP ghost layer
	//if (thisProcData->neighbours[TOP] == MPI_PROC_NULL) {
		//p_setIterationParameters(&endOuter, &endInner, &fixedValue, thisProcData, TOP);
		//for (outerIt = 0; outerIt <= endOuter; outerIt++) {
			//for (innerIt = 0; innerIt <= endInner; innerIt++) {
				//idx = p_computeCellOffset(outerIt,innerIt,fixedValue,thisProcData->xLength, TOP);
				//flagField[idx] = PERIODIC_BOUNDARY;
			//}
		//}
	//}
//}


#include "initLB.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

// int readParameters(int *xlength, double *tau, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
//
//     /* Read values from file given in argv */
//     READ_INT(*argv, *xlength);
//
//     READ_INT(*argv, *timesteps);
//     READ_INT(*argv, *timestepsPerPlotting);
//
//
//     //TODO: (TKS) Adapt this to component case.
//     //if(*tau<=0.5 || *tau>2){
//         //ERROR("Tau is out of stability region (0.5,2.0) (aborting)! \n");
//     //}
//
//   return 0;
// }

void initialiseFields(t_component * c, int *flagField, int xlength){

    /*Setting initial distributions*/
    /*Initializes to equilibrium state with a random density with 1% max from reference density*/

    // current cell index
    int idx, cellIdx;

    // Temporary variables for xlength^2
    int const xlen2 = xlength+2;
    int const xlen2sq = xlen2*xlen2;

    // Temporary variables for z and y offsets
    int zOffset, yzOffset;

    /* initialize collideField and streamField */

    //Intermediate values to calculate feq for a random density.
    double v[3] = {0,0,0};          //Assume no initial velocity

    double rhoVar = 0.01*rhoRef;    //How much initial difference is allowed in density

    int x,y,z;
    srand(5);
    for (int k = 0; k < numComp; k++) {
        for ( z = 1; z <= xlength; ++z) {
            zOffset = z*xlen2sq;
            for ( y = 1; y <= xlength; ++y) {
                yzOffset = y*(xlength+2) + zOffset;
                for ( x = 1; x <= xlength; ++x) {

                    // Compute the base index
                    cellIdx = (yzOffset + x);
                    idx = Q*cellIdx;

                    //Set the initial density to a random offsett to rhoRef
                    double rnd = ((double)rand()/(double)RAND_MAX);
                    c[k].rho[cellIdx] = rhoRef - 0.5*rhoVar + rhoVar*rnd;

                    computeFeqCell(&c[k].rho[cellIdx], v, &c[k].feq[idx]);
                    for (int i = 0; i < Q; ++i) {
                        c[k].collideField[idx+i] = c[k].feq[idx +i];
                        c[k].streamField[idx+i]  = c[k].feq[idx +i];
                    }
                }
            }
        }
    }


    // NOTE: For debug
    // for (z = 1; z <= xlength; z++) {
    //     zOffset = z*xlen2sq;
    //     for (x = 1; x <= xlength; x++) {
    //         idx = Q*(zOffset + 1*(xlength+2) + x);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 1;
    //             c[0].streamField[idx+i] = 1;
    //         }
    //         idx = Q*(zOffset + (xlength)*(xlength+2) + x);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 2;
    //             c[0].streamField[idx+i] = 2;
    //         }
    //     }
    // }
    //
    // for (z = 1; z <= xlength; z++) {
    //     zOffset = z*xlen2sq;
    //     for (y = 1; y <= xlength; y++) {
    //         idx = Q*(zOffset + y*(xlength+2) + 1);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 3;
    //             c[0].streamField[idx+i] = 3;
    //         }
    //         idx = Q*(zOffset + y*(xlength+2) + xlength);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 4;
    //             c[0].streamField[idx+i] = 4;
    //         }
    //     }
    // }
    //
    // for (y = 1; y <= xlength; y++) {
    //     for (x = 1; x <= xlength; x++) {
    //         idx = Q*(xlen2sq + y*(xlength+2) + x);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 5;
    //             c[0].streamField[idx+i] = 5;
    //         }
    //         idx = Q*(xlen2sq*(xlength) + y*(xlength+2) + x);
    //         for (int i = 0; i < Q; i++) {
    //             c[0].collideField[idx+i] = 6;
    //             c[0].streamField[idx+i] = 6;
    //         }
    //     }
    // }
    // y = 1; z = 1;
    // for (x = 1; x <= xlength; x++) {
    //     idx = Q*(z*xlen2sq + y*xlen2 + x);
    //     for (int i = 0; i < Q; i++) {
    //         c[0].collideField[idx + i] = 1;
    //         c[0].streamField[idx + i] = 1;
    //     }
    // }
    //
    // y = xlength; z = xlength;
    // for (x = 1; x <= xlength; x++) {
    //     idx = Q*(z*xlen2sq + y*xlen2 + x);
    //     for (int i = 0; i < Q; i++) {
    //         c[0].collideField[idx + i] = 3;
    //         c[0].streamField[idx + i] = 3;
    //     }
    // }

    // x= 1; y = 1;
    // for (z = 1; z <= xlength; z++) {
    //     idx = Q*(z*xlen2sq + y*xlen2 + x);
    //     for (int i = 0; i < Q; i++) {
    //         c[0].collideField[idx + i] = 1;
    //         c[0].streamField[idx + i] = 1;
    //     }
    // }
    //
    // x = xlength; y = xlength;
    // for (z = 1; z <= xlength; z++) {
    //     idx = Q*(z*xlen2sq + y*xlen2 + x);
    //     for (int i = 0; i < Q; i++) {
    //         c[0].collideField[idx + i] = 3;
    //         c[0].streamField[idx + i] = 3;
    //     }
    // }

    /*Looping over boundary of flagFields*/
    //All points set to zero at memory allocation (using calloc)

    //These are the no-slip walls
    //fixed: z = 0
    for (y = 0; y <= xlength+1; y++) {
        idx = y*(xlength+2);
        for (x = 0; x <= xlength+1; x++) {
            flagField[x+idx] = 1;
        }
    }

    //fixed: x = 0
    //We start at 1 to not include previous cells again from z = 0
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2sq;
        for (y = 0; y <= xlength+1; y++) {
            flagField[zOffset+y*(xlength+2)] = 1;
        }
    }

    //fixed: x = xlength+1
    //We start at 1 to not include previous cells again from z = 0
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2sq + xlength + 1;
        for (y = 0; y <= xlength+1; y++) {
            flagField[zOffset+y*(xlength+2)] = 1;
        }
    }

    //fixed: y = 0
    //from 1:xlength only, to not include cells at upper, lower, left and right edges
    //The edge cells are set in the other loops.
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2sq;
        for (x = 1; x <= xlength; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    //fixed: y = xlength+1
    //same reasoning for index range as in fixed y=0
    for (z = 1; z <= xlength; z++) {
        zOffset = z*xlen2sq + (xlength+1)*(xlength+2);
        for (x = 1; x <= xlength; x++) {
            flagField[zOffset+x] = 1;
        }
    }

    // This is the moving wall. All cells at z=xlength+1 are included (also the edge cells).
    // fixed: z = xlength+1
    zOffset = (xlength+1)*xlen2sq;
    for (y = 0; y <= xlength+1; y++) {
        idx = zOffset + y*(xlength+2);
        for (x = 0; x <= xlength+1; x++) {
            flagField[x+idx] = 2;
        }
    }
}
