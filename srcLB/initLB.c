#include "initLB.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "boundary.h"
#include "computeCellValues.h"
#include <stdio.h>

/*
 * Checks that walls settings are correct. The function is separated into two different cases.
 * The first case checks the correctness of the CAVITY mode, the second case the walls setting
 * for more arbitrary scenarios (scenarios from worksheet 3).
 *
 * The restrictions on the walls are located at comments at the corresponding position in the code.
 */
void p_verifyValidWallSetting(t_boundPara *boundPara, const int mode){

	if(mode == CAVITY){
		//1) Only the 'YZ_TOP wall' is allowed to be the moving wall!
		//2) All the other walls have to be NO_SLIP
		int condBool = 1;

		//1)
		if(boundPara[YZ_TOP].type != MOVING){
			condBool = 0;
		}

		//2)
		for(int b=XY_LEFT; b <= XZ_BACK && condBool; ++b){
			if(b != YZ_TOP && boundPara[b].type != NO_SLIP){
				condBool = 0;
			}
		}

		if(! condBool){
			char msg[150];
			snprintf(msg, 150, "In mode=CAVITY only the YZ_TOP wall is allowed to be a MOVING wall. "
					"All the other walls have to be NO_SLIP");
			ERROR(msg);
		}
	}else{
		int numInflow    = 0;
		int numFreeSlip  = 0;

		// If an INFLOW is present, there must be an OUTFLOW and it has to be on
		// the wall opposite to the INFLOW boundary.
		// We hence also only allow for one inflow in the system.
		//TODO: (TKS) recheck
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

			if(boundPara[i].type == FREE_SLIP || boundPara[i+1].type == FREE_SLIP){
				numFreeSlip++;

			}
		}

		if (numInflow > 1)
			ERROR("Too many inflows");

		if(numFreeSlip > 3)
			ERROR("More than 3 free-slip walls");
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
	boundPara->velocity[0] = x_velocity;
	boundPara->velocity[1] = y_velocity,
	boundPara->velocity[2] = z_velocity;
	boundPara->rhoRef = rhoRef;
	boundPara->rhoIn = rhoIn;

	//Valid setting checking and setting priority of wall (determines order to set)
	int msgSize = 200;
	char msg[msgSize]; // need only in case of errors

	/* Set priorities and check valid settings */
	switch (type) {
		case NO_SLIP:
			boundPara->bPrio = 2;
			break;

		case MOVING:
			boundPara->bPrio = 3;
			break;

		case FREE_SLIP:
			boundPara->bPrio = 1;
			break;

		case INFLOW:
			boundPara->bPrio = 0;

			if(fabs(rhoRef - 1) > densityTol){
				snprintf(msg, msgSize, "Invalid INFLOW wall with rhoRef=%f \n", rhoRef);
				ERROR(msg);
			}
			break;

		case PRESSURE_IN:
			boundPara->bPrio = 0;

			if(fabs(rhoIn - 1) > densityTol){
				snprintf(msg, msgSize, "Invalid PRESSURE_IN with rhoIn=%f \n", rhoIn);
				ERROR(msg);
			}

			break;

		case OUTFLOW:
			boundPara->bPrio = 0;

			if(fabs(rhoRef - 1) > densityTol){
				snprintf(msg, msgSize, "Invalid OUTFLOW wall with rhoRef=%f \n", rhoRef);
				ERROR(msg);
			}

			break;
		default:
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

//TODO: Also include the checking for illegal islands or comment...
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
	//TODO: bring for loop inside checkValidDomain?
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

		if(z_direction > xlength[2] || x_direction > xlength[0] ||
		   z_direction < 0 || x_direction < 0){
            char* str = "";
            sprintf(str, "Domain size = (%d,%d), step size = (%d,%d). Step size cannot be bigger than the domain or negative.",
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
	int MODE = INVALID;
    double Re = INVALID; //u_wall, machNr;
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
    	if((Re == INVALID && *tau == INVALID) || (Re != INVALID && *tau != INVALID)){
        	ERROR("Invalid setting: only provide tau OR Reynolds number (Re). Set the other "
        			"(not used) value to -1 (INVALID).");
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

    for(int b=XY_LEFT; b <= XZ_BACK; b++){
    	printf("\nINFO: reading values for wall type %i \n", b);

    	p_readWall(argv, &boundPara[b], skip);
        skip++;
    }

    p_verifyValidWallSetting(boundPara, MODE);

    /* S T A R T  -- ONLY FOR CAVITY, this part supports the old CAVITY case from WS2
     * Also it computes tau from the given Reynolds number*/

    //CAVITY supports to hand either Reynolds number (and compute tau), or handle tau directly
    if(Re != INVALID && MODE == CAVITY){ //true, then compute tau according to Re
    	//The restriction of cubic domain is required in CAVITY such that the
    	//computation of tau is valid.
    	if(xlength[0] != xlength[1] && xlength[0] != xlength[2]){
    		ERROR("Only cubic domain is supported for the CAVITY mode. \n");
    	}

    	//From p_verifyValidWallSetting it is only possible that the YZ_TOP wall is MOVING
    	assert(boundPara[YZ_TOP].type == MOVING);

    	double u_wall = sqrt(boundPara[YZ_TOP].velocity[0]*boundPara[YZ_TOP].velocity[0]
			+ boundPara[YZ_TOP].velocity[1]*boundPara[YZ_TOP].velocity[1]+
			boundPara[YZ_TOP].velocity[2]*boundPara[YZ_TOP].velocity[2]);

    	*tau = u_wall*(xlength[0])/(C_S*C_S*Re)+0.5;
    	printf("\nINFO: Calculated tau = %f \n", *tau);

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

/* Sets the values for according to 'wallPos'. Because the offsets have to be computed
 * differently for each 'wallPos' the loop indices are first selected.
 */
void p_setWall(const int * const xlength, t_flagField *flagField,
				const int wallType, const int wallPos) {
	int inner, outer, flagIndex;
	int index[2] = {0};
	int point[3] = {0};

	//select looping indices:
	if (wallPos == XY_LEFT) {
		index[0] = 0;
		index[1] = 1;
		point[2] = 0;
	} else if (wallPos == XY_RIGHT) {
		index[0] = 0;
		index[1] = 1;
		point[2] = xlength[2]+1;
	} else if (wallPos == XZ_BACK) {
		index[0] = 0;
		index[1] = 2;
		point[1] = 0;
	} else if (wallPos == XZ_FRONT) {
		index[0] = 0;
		index[1] = 2;
		point[1] = xlength[1]+1;
	} else if (wallPos == YZ_BOTTOM) {
		index[0] = 1;
		index[1] = 2;
		point[0] = 0;
	} else if (wallPos == YZ_TOP) {
		index[0] = 1;
		index[1] = 2;
		point[0] = xlength[0]+1;
	} else {
		ERROR("This should not happen");
	}

	//loop over wall and set accordingly
	for (outer = 0; outer < xlength[index[1]]+2; outer++) {
		for (inner = 0; inner < xlength[index[0]]+2; inner++) {
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

    //Set boundary cells (ghost layer first)
    //Set the walls in order of the wall priorities. The priorities are set in p_readWall
	for (int i = 0; i <= 3; i++) {
		for (int j = 0; j < NUM_WALLS; j++) {
			if (boundPara[j].bPrio == i) {
				p_setWall(xlength, flagField, boundPara[j].type, j);
			}
		}
	}

    //NOTE: only domain (ghost layer were set previously)//
    for(z = 1; z <= xlength[2]; ++z){
    	for (y = 1; y <= xlength[1]; ++y) {
			for (x = 1; x <= xlength[0]; ++x) {
				p_computeIndexXYZ(x,y,z,xlength,&idx);
				int type_domain = pgmMatrix[z][x];

				if(type_domain == 0){
					flagField[idx].type = FLUID;
					flagField[idx].position = FLUID; //the position is needed for boundaries primarily
				}else if(type_domain == 1){          //therefore here it is just set to corresponding type
					flagField[idx].type = OBSTACLE;
					flagField[idx].position = OBSTACLE;
				}else{
					ERROR("Description of scenario in pgm file should only consist of logical values. \n");
				}
			}
		}
    }

    /* initialize collideField and streamField*/
    //f_i(x,0) = f^eq(1,0,0) = w_i
    for ( z = 0; z < zlen2; ++z) {
        for ( y = 0; y < ylen2; ++y) {
            for ( x = 0; x < xlen2; ++x) {
				// Compute the cell index
                p_computeIndexXYZ(x,y,z,xlength,&idx);
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
