#include "initLB.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "boundary.h"
#include "computeCellValues.h"
#include <stdio.h>



/*
 * Checks that walls settings are correct. The function is separated into two different cases.
 * The first case checks the correctness of the CAVITY mode, the second case the walls setting
 * for more arbitrary scenarios (scenarios from worsheet 3).
 *
 * The restrictions on the walls are located at comments at the corresponding position in the code.
 */
void p_verifyValidWallSetting(t_boundPara *boundPara, const int mode){

	if(mode == CAVITY){
		//If there is a MOVING wall present we only support the cavity scenario, i.e. there should
		//be only one MOVING wall and 5 NO_SLIP walls.
		int foundMovWall = 0;
		int foundNoSlipWall = 0;
		for(int b=XY_LEFT; b <= XY_LEFT; ++b){
			if(boundPara->type == MOVING){
				foundMovWall++;
			}
			if(boundPara->type == NO_SLIP){
				foundNoSlipWall++;
			}
		}

		if(foundMovWall > 1 || (foundMovWall == 1 && foundNoSlipWall != 5)){
			char msg[150];
			snprintf(msg, 150, "When a wall is set to MOVING we only support the CAVITY scenario. That is "
					"one MOVING wall and 5 NO_SLIP walls.");
			ERROR(msg);
		}// else if(foundMovWall == 1 && foundNoSlipWall == 5){

	}else{
		int numInflow    = 0;
		int numFreeSlip  = 0;

		// If an INFLOW is present, there must be an OUTFLOW and it has to be on
		// the wall opposite to the INFLOW boundary.
		// We hence also only allow for one inflow in the system.
		for(int i = XY_LEFT; i <= XZ_BACK; i=i+2) {
			if(boundPara[i].type == MOVING || boundPara[i+1].type == MOVING){
				ERROR("MOVING walls are only supported in mode CAVITY!!");
			}

			if(boundPara[i].type == INFLOW){
				if(boundPara[i+1].type != OUTFLOW){
					ERROR("OUTFLOW is not opposite INFLOW");
				}
				numInflow++;
			}

			if(boundPara[i].type == OUTFLOW){
				if(boundPara[i+1].type != INFLOW){
					ERROR("OUTFLOW is not opposite INFLOW");
				}
				numInflow++;
			}

			if(boundPara[i].type == FREE_SLIP){
				numFreeSlip++;
				if (numFreeSlip == 2 && (boundPara[i+1].type != FREE_SLIP)) {
					ERROR("FREE_SLIP walls are not opposites");
				}
			}
			if(boundPara[i+1].type == FREE_SLIP){
				numFreeSlip++;
				if (numFreeSlip == 2 && (boundPara[i].type != FREE_SLIP)) {
					ERROR("FREE_SLIP walls are not opposites");
				}
			}
		}

		if (numInflow > 1)
			ERROR("Too many inflows");

		if(numFreeSlip > 3)
			ERROR("Too many free-slip walls");
	}
}


/*
 * Reads a single wall from the parameter file that is provided. Additionally, some checks
 * for correctness of the user-settings are carried out.
 */
void p_readWall(char *argv[], t_boundPara *boundPara, const int skip){

	int type;
	double x_velocity, y_velocity, z_velocity;
	double rhoRef, rhoIn;

    READ_INT(*argv, type, skip);

	READ_DOUBLE(*argv, x_velocity, skip);
	READ_DOUBLE(*argv, y_velocity, skip);
    READ_DOUBLE(*argv, z_velocity, skip);

	READ_DOUBLE(*argv, rhoRef, skip);
	READ_DOUBLE(*argv, rhoIn, skip);

	boundPara->type = type;
	boundPara->wallVelocity[0] = x_velocity;
	boundPara->wallVelocity[1] = y_velocity,
	boundPara->wallVelocity[2] = z_velocity;
	boundPara->rhoRef = rhoRef;
	boundPara->rhoIn = rhoIn;

	//Valid setting checking:
	int msgSize = 200;
	char msg[msgSize]; // only needed in case of errors


	/* TODO: (DL) discuss tolerances that are allowed */
	switch(type){
	case INFLOW:
		if(fabs(rhoRef - 1) > 0.3){
			snprintf(msg, msgSize, "Invalid INFLOW wall with rhoRef=%f \n", rhoRef);
			ERROR(msg);
		}
		break;
	case OUTFLOW:
		if(fabs(rhoRef - 1) > 0.3){
			snprintf(msg, msgSize, "Invalid OUTFLOW wall with rhoRef=%f \n", rhoRef);
			ERROR(msg);
		}
		break;
	case PRESSURE_IN:
		if(fabs(rhoIn - 1) > 0.2){
			snprintf(msg, msgSize, "Invalid PRESSURE_IN with rhoIn=%f \n", rhoIn);
			ERROR(msg);
		}
		break;
	}

	if(boundPara->type>=NO_SLIP && boundPara->type<=MOVING){ //NO_SLIP OR MOVING
		boundPara->idxStartEnd[0] = 0;
		boundPara->idxStartEnd[1] = 1;

	}else if(boundPara->type>=FREE_SLIP && boundPara->type<=PRESSURE_IN){ //FREE_SLIP, INFLOW, OUTFLOW, PRESSURE_IN
		boundPara->idxStartEnd[0] = 1;
		boundPara->idxStartEnd[1] = 0;
	}else{
		snprintf(msg, msgSize, "Type %i, is not valid. See E_CELL_TYPE in LBDefinitions.h for valid"
				" types.", type);
		ERROR(msg);
	}
}

/*
* Checks whether the surroundings of the current cell is legal.
* Geometry that is not allowed:
*  ##       #
*    ##     #   ("thin wall")
*           #
*           #
* Returns an boolean value to indicate whether the surrounding of this position is valid.
* The function only checks custom scenarios and is thus only called in function
* p_readCustomPgmMatrix
*/
int p_checkValidDomain(int x, int z, const int* const xlength, int **pgmMatrix){

	int right = pgmMatrix[z+1][x];
    int left  = pgmMatrix[z-1][x];
    int up    = pgmMatrix[z][x+1];
    int down  = pgmMatrix[z][x-1];
    int is_valid = 1;

    if( !(right || up)){                    // If right and up are FLUID
        if(pgmMatrix[z+1][x+1] == 1) { 		// If right corner is an OBSTACLE
            is_valid = 0;
        }
    }

    if( !(right || down)){                  // If right and down are FLUID
        if(pgmMatrix[z+1][x-1] == 1) { 		// If right corner is an OBSTACLE
            is_valid = 0;
        }
    }

    if( !(left || up)){                     // If left and up are FLUID
        if(pgmMatrix[z-1][x+1] == 1) { 		// If left corner is an OBSTACLE
            is_valid = 0;
        }
    }

    if( !(left || down)){                   // If left and down are FLUID
        if(pgmMatrix[z-1][x-1] == 1) { 		// If left corner is an OBSTACLE
            is_valid = 0;
        }
    }
    return(is_valid);
}

/*Global variable makes it cleaner. A domain is described as in a pgm format
 * even when no pgm file is provided. The two functions 'readParameters' and
 * 'initialiseFields' use it:
 * 'readParameters' - allocates and checks for validity of the geometry
 * 'initialiseFields' - uses the matrix to set the domain and deallocates the matrix
*/
static int **pgmMatrix = NULL;

/*
 * Allocation of the pgmMatrix in cases where the domain is created according to parameters provided,
 * i.e. no extern pgm-file.
 * The function sets the values for the boundary layer (the outer values are 1) and the domain
 * values are 0. (Only logical values are used for the pgm matrix)
 */
void p_initPgmMatrix(const int * const xlength){

	pgmMatrix = imatrix(0, xlength[2]+1, 0, xlength[0]+1); //XZ-plane where Z is on the "x-axis"

	// set all values to 1 in the matrix
	init_imatrix(pgmMatrix, 0, xlength[2]+1, 0, xlength[0]+1, 1);

	// set all domain values to 0 (leaving the boundary to 1)
	init_imatrix(pgmMatrix, 1, xlength[2], 1, xlength[0], 0);
}

/*
 * Allocation of pgmMatrix for custom scenarios. If the domain sizes (*xlength) of the provided
 * parameter file and the pgm matrix do not correspond an error message will show up.
 *
 * Furthermore, the domain is checked if it is valid (see also function p_checkValidSurrounding)
 */
void p_readCustomPgmMatrix(const int * const xlength, char *filename){
	char file_path[MAX_LINE_LENGTH+10]; //+10 for 'scenarios/'
	snprintf(file_path, MAX_LINE_LENGTH+10, "scenarios/%s", filename);
	int zsizePgm, xsizePgm;
	pgmMatrix = read_pgm(file_path, &zsizePgm, &xsizePgm);

	if(zsizePgm != xlength[2] || xsizePgm != xlength[0]){
		free_imatrix(pgmMatrix,0,zsizePgm+1,0,xsizePgm+1);
		char msg[150];
		snprintf(msg, 150, "The size of the pgm matrix is [z_length=%i x xlength=%i]. The parameter"
				"read in the setting file are: [z_length=%i x xlength=%i] \n", zsizePgm, xsizePgm, xlength[2], xlength[0]);
		ERROR(msg);
	}

	/*Check for illegal geometries*/
    for(int z = 1; z <= xlength[2]; ++z){
        for (int x = 1; x <= xlength[0]; ++x) {
            if(pgmMatrix[z][x] == 1){
                if(!p_checkValidDomain(x, z, xlength, pgmMatrix)){
		            free_imatrix(pgmMatrix,0,zsizePgm+1,0,xsizePgm+1);
                    char error[80];
                    snprintf(error, 80, "Invalid surroundings at x = %d, z = %d)", x,z);
                    ERROR(error);
                }
			}
		}
    }
}

/* TODO: (DL) delete when finializing code, but leave it as long there may be debugging
 * going on :-) */
void print_matrix(const int * const xlength){
	 for (int x = 0; x <= xlength[0]+1; x++) {
       for (int z = 0; z <= xlength[2]+1; z++) {
           printf("%d ",pgmMatrix[z][xlength[0]+1-x]);
       }
       printf("\n");
   }
}

/*
 * Handles the creation (and allocation) of the pgmMatrix. The domain is created depending
 * on the current mode.
 */
void p_generatePgmDomain(char *argv[], const int * const xlength, const int MODE){

	//cannot be allocated in "case" -statement (compiler error)
	//it's only needed in the arbitrary case
	char problem[MAX_LINE_LENGTH];

	switch(MODE){
	case SHEAR_FLOW:
		p_initPgmMatrix(xlength);
		break;
	case CAVITY:
		p_initPgmMatrix(xlength);
		break;
	case STEP_FLOW:
		p_initPgmMatrix(xlength);
		//insert step:
		int z_direction, x_direction;
		READ_INT(*argv, z_direction,0); //size of step
		READ_INT(*argv, x_direction,0);

		if(z_direction > xlength[2] || x_direction > xlength[0]){
            char* str = "";
            sprintf(str, "Domain size = (%d,%d), step size = (%d,%d). Step size cannot be bigger than the domain.",
                          xlength[0], xlength[2],x_direction,z_direction);
			ERROR(str);
		}

		init_imatrix(pgmMatrix, 1, z_direction, 1, x_direction, 1);
		break;
	case ARBITRARY:
		READ_STRING(*argv, problem, 0);
		p_readCustomPgmMatrix(xlength, problem);
		break;
	default:
		ERROR("SETTING-MODE is invalid!!");
		break;
	}
}


int readParameters(int *xlength, double *tau, t_boundPara *boundPara, int *timesteps,
		int *timestepsPerPlotting, char *problem,
		int argc, char *argv[]){

	if(argc != 2){
		char msg[200];
		snprintf(msg, 200, "There are %i arguments provided to the simulation. Only 1 argument "
				"providing \nthe path+filename (of relative to working directory) of the scenario "
				"is accepted!", argc-1);
		ERROR(msg);
	}

	printf("INFO: read simulation parameter: \n");
	int MODE;
    double Re; //u_wall, machNr;
    int x_length, y_length, z_length; //Temporary for reading
    int skip = 0; // How many of the same named variable should be skipped when
                  // calling READ_<TYPE>.

    READ_INT(*argv, MODE, 0);
    if(MODE < 0 || MODE >= NUM_MODES ){
    	char msg[150];
    	snprintf(msg, 150, "SETTING-MODE=%i is invalid! Look at enum SETTING_MODE in LBDefintions.h"
    			"for valid modes. \n", MODE);
    	ERROR(msg);
    }

    READ_DOUBLE(*argv, *tau, 0);

    if(MODE == CAVITY){
        READ_DOUBLE(*argv, Re, 0);
    	if((Re == -1 && *tau == -1) || (Re != -1 && *tau != -1)){
        	ERROR("Invalid setting: only provide tau OR Reynolds number (Re). Set the other "
        			"(not used) value to -1.");
        }
    }

    READ_INT(*argv, *timesteps, 0);
    READ_INT(*argv, *timestepsPerPlotting, 0);

    //pgm-file OR simply the name that describes the scenario
    //pgm file has to be located in /scenarios
    READ_STRING(*argv, problem, 0);

    //Domain size
    READ_INT(*argv, x_length, 0);
    READ_INT(*argv, y_length, 0);
    READ_INT(*argv, z_length, 0);
    xlength[0] = x_length;
    xlength[1] = y_length;
    xlength[2] = z_length;

    //TODO: (TKS) Make sure the boundaries are red in the correct order.
    for(int b=XY_LEFT; b <= XZ_BACK; b++){
    	printf("\nINFO: reading values for wall type %i \n", b);

    	p_readWall(argv, &boundPara[b], skip);
        skip++;
    }

    p_verifyValidWallSetting(boundPara, MODE);

    /* S T A R T  -- ONLY FOR CAVITY, this part supports the old CAVITY case from WS2
     * Also it computes tau from the given Reynolds number
     */

    //CAVITY supports to hand either Reynolds number (and compute tau), or handle tau directly
    if(Re != -1.0 && MODE == CAVITY){ //true, then compute tau according to Re
    	//The restriction of cubic domain is required in CAVITY such that the
    	//computation of tau is valid.
    	if(xlength[0] != xlength[1] && xlength[0] != xlength[2]){
    		ERROR("Only quadratic domain is supported for the CAVITY mode. \n");
    	}
    	double u_wall = -1;
    	for(int b=XY_LEFT; b <= XZ_BACK; ++b){
    	    if(boundPara[b].type == MOVING){
    			u_wall = sqrt(boundPara[b].wallVelocity[0]*boundPara[b].wallVelocity[0]
					+ boundPara[b].wallVelocity[1]*boundPara[b].wallVelocity[1]+
					boundPara[b].wallVelocity[2]*boundPara[b].wallVelocity[2]);
    			*tau = u_wall*(xlength[0])/(C_S*C_S*Re)+0.5;
    			printf("\nINFO: Calculated tau = %f \n", *tau);

    			break; //CAVITY case, in case there are multiple MOVING walls an
    				   //error is thrown in valid surroundings
    		}
    	}
    	*tau = u_wall*(*xlength)/(C_S*C_S*Re)+0.5;

    	if(u_wall == -1) ERROR("Couldn't find a MOVING wall.");

    	double machNr  = u_wall/C_S;

    	printf("\nINFO: Wall speed = %f \n", u_wall);
    	printf("\nINFO: Mach number = %f \n\n", machNr);

    	/* valid settings check*/
    	if(u_wall >= C_S){
    		ERROR("Wall speed is supersonic (aborting)! \n");
    	}

    	/*We allow user defined mach number tolerance for Ma << 1 (default = 0.1)
    	      To change please look at LBDefinitions.h*/
    	if(machNr > machNrTol){
    		char buffer[80];
    		snprintf(buffer, 80, "Mach number is larger than %f (aborting)! \n",machNrTol);
    		ERROR(buffer);
    	}
    }
    /* E N D  -- CAVITY PART*/


    if(*tau<=0.5 || *tau>2){
    	char msg[80];
    	snprintf(msg, 80, "Tau (value=%f) is out of stability region (0.5,2.0). ! \n", *tau);
    	ERROR(msg);
    }

    p_generatePgmDomain(argv, xlength, MODE);

    return 0;
}
//TODO: (TKS) REMOVE AFTER TESTING
void print_flagfield_slice(int* field, const int * const xlength){
    int idx;
    int y = xlength[1]/2;

    for (int x = 0; x <= xlength[0]+1; x++) {
        for (int z = 0; z <= xlength[2]+1; z++) {
            p_computeIndexXYZ(x,y,z,xlength,&idx);
            printf("%d ",field[idx]);
        }
        printf("\n");
    }
}

void p_setWall(const int * const xlength, const int fixed, const int outerIdx, const int innerIdx, t_flagField *flagField,
				const int wallType, const int wallPos) {
	int inner, outer, flagIndex;
	int index[2];
	int point[3];
	if (wallPos == XY_LEFT || wallPos == XY_RIGHT) {
		index[0] = 0;
		index[1] = 1;
		point[2] = fixed;
	} else if (wallPos == XZ_BACK || wallPos == XZ_FRONT) {
		index[0] = 0;
		index[1] = 2;
		point[1] = fixed;
	} else if (wallPos == YZ_BOTTOM || wallPos == YZ_TOP) {
		index[0] = 1;
		index[1] = 2;
		point[0] = fixed;
	} else {
		ERROR("This should not happen");
	}
	for (outer = 0; outer < xlength[outerIdx]+2; outer++) {
		for (inner = 0; inner < xlength[innerIdx]+2; inner++) {
			point[index[0]] = inner;
			point[index[1]] = outer;
			p_computeIndex(point, xlength,&flagIndex);
			flagField[flagIndex].type = wallType;
			flagField[flagIndex].position = wallPos;
		}
	}
}

void initialiseFields(double *collideField, double *streamField, t_flagField *flagField,
		int *xlength, t_boundPara *boundPara, char *problem){

	// Loop variables
	int x,y,z;

    // Temporary variables for xlength^2
	int const xlen2		= xlength[0]+2;
	int const ylen2		= xlength[1]+2;
	int const zlen2		= xlength[2]+2;
	int idx;

    int type;
    int start, end;

    //Do the domain enclosing boundaries (ghost layers) first and set type accordingly.

	// y = 0
    type = boundPara[XZ_BACK].type;
    start = boundPara[XZ_BACK].idxStartEnd[0];
    end = boundPara[XZ_BACK].idxStartEnd[1];
    for(z = start; z <= xlength[2] + end; ++z){
    	for(x = start; x <= xlength[0] + end; ++x){
			p_computeIndexXYZ(x,0,z,xlength,&idx);
    		flagField[idx].type = type;
			flagField[idx].position = XZ_BACK;
    	}
    }

	// y = xlength[1]+1
    type = boundPara[XZ_FRONT].type;
    start = boundPara[XZ_FRONT].idxStartEnd[0];
    end = boundPara[XZ_FRONT].idxStartEnd[1];
    for(z = start; z <= xlength[2]+end; ++z){
    	for(x = start; x <= xlength[0]+end; ++x){
			p_computeIndexXYZ(x,xlength[1]+1,z,xlength,&idx);
    		flagField[idx].type = type;
			flagField[idx].position = XZ_FRONT;
    	}
    }

	// z = 0
    type = boundPara[XY_LEFT].type;
    start = boundPara[XY_LEFT].idxStartEnd[0];
    end = boundPara[XY_LEFT].idxStartEnd[1];
    for(y = start; y <= xlength[1]+end; ++y){
    	for(x = start; x <= xlength[0]+end; ++x){
			p_computeIndexXYZ(x,y,0,xlength,&idx);
    		flagField[idx].type = type;
			flagField[idx].position = XY_LEFT;
    	}
    }

	// z = xlength[2]+1
    type = boundPara[XY_RIGHT].type;
    start = boundPara[XY_RIGHT].idxStartEnd[0];
    end = boundPara[XY_RIGHT].idxStartEnd[1];
    for(y = start; y <= xlength[1]+end; ++y){
    	for(x = start; x <= xlength[0]+end; ++x){
			p_computeIndexXYZ(x,y,xlength[2]+1,xlength,&idx);
    		flagField[idx].type = type;
			flagField[idx].position = XY_RIGHT;
    	}
    }

	// x = 0
    type = boundPara[YZ_BOTTOM].type;
    start = boundPara[YZ_BOTTOM].idxStartEnd[0];
    end = boundPara[YZ_BOTTOM].idxStartEnd[1];
    for(z = start; z <= xlength[2]+end; ++z){
    	for(y = start; y <= xlength[1]+end; ++y){
			p_computeIndexXYZ(0,y,z,xlength,&idx);
			flagField[idx].type = type;
			flagField[idx].position = YZ_BOTTOM;
    	}
    }

	// x = xlength[0]+1
    type = boundPara[YZ_TOP].type;
    start = boundPara[YZ_TOP].idxStartEnd[0];
    end = boundPara[YZ_TOP].idxStartEnd[1];
    for(z = start; z <= xlength[2]+end; ++z){
    	for(y = start; y <= xlength[1]+end; ++y){
			p_computeIndexXYZ(xlength[0]+1,y,z,xlength,&idx);
    		flagField[idx].type = type;
			flagField[idx].position = YZ_TOP;
    	}
    }

    /* TODO: (DL) only for testing START */
//    for(y = 0; y <= xlength[1]+1; ++y){
//		offset1 = y*xlen2;
//		offset2 = offset1 + (xlength[2]+1)*xlen2*ylen2;
//    	for(x = 0; x <= xlength[0]+1; ++x){
//
////    		printf("%i ", flagField[offset1 + x].type);
//
//    		if(flagField[offset1 + x].type == -1){ // z = 0
//    			flagField[offset1 + x].type = PRESSURE_IN;
//    			flagField[offset1 + x].position = XY_LEFT;
//    		}
//    	}
////    	printf(" \n");
//    }
//
//    type = boundPara[XY_RIGHT].type;
//    for(y = 0; y <= xlength[1]+1; ++y){
//		offset1 = y*xlen2;
//		offset2 = offset1 + (xlength[2]+1)*xlen2*ylen2;
//    	for(x = 0; x <= xlength[0]+1; ++x){
//
//    		if(flagField[offset2 + x].type == -1){ // z = 0
//    			flagField[offset2 + x].type = OUTFLOW;
//    			flagField[offset2 + x].position = XY_RIGHT;
//    		}
//    	}
//    }
      /* TODO: (DL) only for testing END */

    //NOTE: only domain (ghost layer were set previously)//
    for(z = 1; z <= xlength[2]; ++z){
    	for (y = 1; y <= xlength[1]; ++y) {
			for (x = 1; x <= xlength[0]; ++x) {
				p_computeIndexXYZ(x,y,z,xlength,&idx);
				int type_domain = pgmMatrix[z][x];

				if(type_domain == 0){
					// TODO: (REMOVE?)This is actually not required as long as FLUID=0
					flagField[idx].type = FLUID;
					flagField[idx].position = FLUID;
				}else if(type_domain == 1){
					flagField[idx].type = OBSTACLE;
					flagField[idx].position = OBSTACLE;
				}else{
					ERROR("Description of scenario in pgm file should only consist of logical values. \n");
				}
			}
		}
    }

    //TODO: (TKS) remove when finished testing
    //print_flagfield_slice(flagField, xlength);

    /*Setting initial distributions, LATTICEWEIGHTS and INFLOW conditions */
    //f_i(x,0) = f^eq(1,0,0) = w_i

    /* initialize collideField and streamField and inflow BC if any*/
    for ( z = 0; z < zlen2; ++z) {
        for ( y = 0; y < ylen2; ++y) {
            for ( x = 0; x < xlen2; ++x) {
				p_computeIndexXYZ(x,y,z,xlength,&idx);
				// Compute the base index
                idx = Q*idx;

                //Set initial condition
                for (int i = 0; i < Q; ++i) {
                    collideField[idx+i] = LATTICEWEIGHTS[i];
                    streamField[idx+i]  = LATTICEWEIGHTS[i];
                }
            }
        }
    }

    free_imatrix(pgmMatrix,0,xlength[2]+1,0,xlength[0]+1);
}
