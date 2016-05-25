#include "initLB.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "boundary.h"
#include "computeCellValues.h"
#include <stdio.h>

// Function for setting the inflow once and for all
void p_handleInflow(int x, int y, int z, int *xlength, t_boundPara *boundPara,
                     double *collideField, const int currentCellIndex){
	int i, flag;
	double feq[19];
	if(x == 0)
		flag = YZ_BOTTOM;
	else if(x == xlength[0]+1)
		flag = YZ_TOP;
	else if(y == 0)
		flag = XZ_BACK;
	else if(y == xlength[1]+1)
		flag = XZ_FRONT;
	else if(z == 0)
		flag = XY_LEFT;
	else if(z == xlength[2]+1)
		flag = XY_RIGHT;
	else{
		flag = -1;
		ERROR("Inflow not a boundary (Aborting)");
	}
	if (flag != -1) {
		computeFeq(&boundPara[flag].rhoRef, boundPara[flag].wallVelocity,feq);
		for(i = 0; i < Q; i++){
			for(int i = 0; i < Q; i++){
				collideField[currentCellIndex + i] = feq[i];
			}
		}
	}
 }

/*
TODO: (DL) SUGGESTION FOR INDEXING
A "standard" scenario is that the pipe goes from left to right. So the inflow should really only happen where there is no boundary starting.

This would be a bad setting (in 2D): I - inflow, F - fluid, B -boundary

I F F
I F F
I B B

In the examples the YZ plane is set to no-slip - meaning a propper boundary - so it should cover the entire plane
By that we can reconstruct our cavity also in the same way we had in the last WS because we need the top all MOVING.

The other question is how we set the overlapping cells of the XZ plane (in the examples often free-slip)

I guess the same argument holds:

S-Free slip
I F F
I F F
I S S

this looks kind of bad and I have the feeling this could cause troubles...

So vote that the XZ plane should get the cells that overlap with the XY plane.
By that the (normally) inflow cells are only streaming directly towards fluid cells (and not diagonal from an edge)

And I think the more important thing is, that our examples look nice, and not that _every_ possibility
looks nice ;-)

By these settings the examples make more sense!
 */

void p_verifyValidWallSetting(t_boundPara *boundPara){

    int num_inflow    = 0;
    int num_free_slip = 0;

    // If an INFLOW is present, there must be an OUTFLOW and it has to be on
    // the wall oppopsite to the INFLOW.
    // We hence also only allow for one inflow on the boundary.

    for (int i = XY_LEFT; i <= XZ_BACK; i=i+2) {
        if(boundPara[i].type == INFLOW){
            if(boundPara[i+1].type != OUTFLOW){
                ERROR("OUTFLOW is not opposite INFLOW");
            }
            num_inflow++;
        }

        if(boundPara[i].type == OUTFLOW){
            if(boundPara[i+1].type != INFLOW){
                ERROR("OUTFLOW is not opposite INFLOW");
            }
            num_inflow++;
        }

        if(boundPara[i].type == FREE_SLIP){
            num_free_slip++;
            if (num_free_slip == 2 && (boundPara[i+1].type != FREE_SLIP)) {
                    ERROR("FREE_SLIP walls are not opposites");
            }
        }
        if(boundPara[i+1].type == FREE_SLIP){
            num_free_slip++;
            if (num_free_slip == 2 && (boundPara[i].type != FREE_SLIP)) {
                    ERROR("FREE_SLIP walls are not opposites");
            }
        }
    }

    if (num_inflow > 1)
        ERROR("Too many inflows");

    if(num_free_slip > 3)
        ERROR("Too many free-slip walls");
}



void p_readWall(char *argv[], t_boundPara *boundPara, const int skip){
	int type;
	double x_velocity, y_velocity, z_velocity;
	double rhoRef, rhoIn; /* TODO: (DL) rename rhoIn -> deltaRho, also in templates*/

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


	/* TODO: (DL) after testing this can be much simplified (only in two ifs), but for testing
	 * this version is a bit more handy...
	 */
	switch (boundPara->type) {
	case NO_SLIP:
		boundPara->idxStartEnd[0] = 0;
		boundPara->idxStartEnd[1] = 1;
		break;
	case MOVING:
		boundPara->idxStartEnd[0] = 0;
		boundPara->idxStartEnd[1] = 1;
		break;
	case FREE_SLIP:
		boundPara->idxStartEnd[0] = 1;
		boundPara->idxStartEnd[1] = 0;
		break;
	case INFLOW:
		boundPara->idxStartEnd[0] = 1;
		boundPara->idxStartEnd[1] = 0;
		break;
	case OUTFLOW:
		boundPara->idxStartEnd[0] = 1;
		boundPara->idxStartEnd[1] = 0;
		break;
	case PRESSURE_IN:
		boundPara->idxStartEnd[0] = 1;
		boundPara->idxStartEnd[1] = 0;
		break;
	default:
		//TODO: (DL) msg
		ERROR("TODO: error msg");
	}
}

//Checks whether the surroundings of the current cell is legal.
//Examples of geometries not allowed:
//  ##       #
//    ##     #
//            #
//            #

int valid_sorroundings(int x, int z, const int* const xlength, int **pgmMatrix){

    int right = pgmMatrix[z+1][x];
    int left  = pgmMatrix[z-1][x];
    int up    = pgmMatrix[z][x+1];
    int down  = pgmMatrix[z][x-1];
    int is_valid = 1;

    if( !(right && up)){                    // If right and up are FLUID
        if(pgmMatrix[x+1][z+1] == 1) 		// If right corner is an OBSTACLE
            is_valid = 0;
    }

    if( !(right && down)){                  // If right and down are FLUID
        if(pgmMatrix[x-1][z+1] == 1) 		// If right corner is an OBSTACLE
            is_valid = 0;
    }

    if( !(left && up)){                     // If left and up are FLUID
        if(pgmMatrix[x+1][z-1] == 1)        // If right corner is an OBSTACLE
            is_valid = 0;
    }

    if( !(left && down)){                   // If left and down are FLUID
        if(pgmMatrix[x-1][z-1] == 1)        // If right corner is an OBSTACLE
            is_valid = 0;
    }

    return(is_valid);

}

//Global variable makes it cleaner. A domain is described as in a pgm format
//even when no pgm file is provided. The two functions 'readParameters' and
//'initialiseFields' use it:
//'readParameters' - allocates and checks for validity of the geometry
//'initialiseFields' - uses the matrix to set the domain and deallocates the matrix
static int **pgmMatrix = NULL;

void print_matrix(const int * const xlength);

void init_pgmMatrix(const int * const xlength){

	pgmMatrix = imatrix(0, xlength[2]+1, 0, xlength[0]+1); //XZ-plane where Z is on the "x-axis"

	// set all values to 1 in the matrix
	init_imatrix(pgmMatrix, 0, xlength[2]+1, 0, xlength[0]+1, 1);

	// set all domain values to 0 (leaving the boundary to 1)
	//TODO: (DL) when I comment out this line in the, the freeing does not encounter any errors.
	init_imatrix(pgmMatrix, 1, xlength[2], 1, xlength[0], 0);
}

void read_customPgmMatrix(const int * const xlength, char *filename){
	char file_path[MAX_LINE_LENGTH+10]; //+10 for 'scenarios/'
	snprintf(file_path, MAX_LINE_LENGTH+10, "scenarios/%s", filename);
	int zsizePgm, xsizePgm;
	pgmMatrix = read_pgm(file_path, &zsizePgm, &xsizePgm);

	if(zsizePgm != xlength[2] && xsizePgm != xlength[0]){
		free_imatrix(pgmMatrix,0,zsizePgm+1,0,xsizePgm+1);
		ERROR("TODO: msg");
	}

    /*Check for illegal geometries*/
    for(int z = 1; z <= xlength[2]; ++z){
        for (int x = 1; x <= xlength[0]; ++x) {
            if(pgmMatrix[z][x] == 1){
                if(!valid_sorroundings(x, z, xlength, pgmMatrix)){
		            free_imatrix(pgmMatrix,0,zsizePgm+1,0,xsizePgm+1);
                    char* error = "";
                    sprintf(error, "Invalid surroundings at x = %d, z = %d)", x,z);
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

void p_generatePgmDomain(char *argv[], const int * const xlength, const int MODE){

	//cannot be allocated in "case" -statement (compiler error)
	//it's only needed in the arbitrary case
	char problem[MAX_LINE_LENGTH];

	switch(MODE){

	case SHEAR_FLOW:
		init_pgmMatrix(xlength);
		break;
	case CAVITY:
		init_pgmMatrix(xlength);
		break;
	case STEP_FLOW:
		init_pgmMatrix(xlength);
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
		read_customPgmMatrix(xlength, problem);
		break;
	default:
		ERROR("SETTING-MODE is invalid!!");
		break;
	}

}

// TODO: (TKS) The printed statements from the helper function do no longer
//             make sense after adding skip
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

	int MODE;
    double Re; //u_wall, machNr;
    int x_length, y_length, z_length; //Temporary for reading
    int skip = 0; // How many of the same named variable should be skipped when
              // calling READ_<TYPE>.

    READ_INT(*argv, MODE, 0);
    if(MODE < 0 || MODE >= NUM_MODES ){
    	char msg[40];
    	snprintf(msg, 40, "SETTING-MODE=%i is invalid! \n", MODE);
    	ERROR(msg);
    }

    READ_DOUBLE(*argv, Re, 0);
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
    	p_readWall(argv, &boundPara[b], skip);
        skip++;
    }
    p_verifyValidWallSetting(boundPara);

    p_generatePgmDomain(argv, xlength, MODE);

    /*TODO: (DL) the characteristic velocity, characteristic length, mach number are
     * no longer valid (using values valid for cavity)
     *
     * I suggest we create a new function where we check for valid settings only!
     * For now all checks are comment out.
     */

    /*Calculates tau from the Reynolds number*/

//    u_wall  = sqrt(xvelocity*xvelocity + yvelocity*yvelocity+zvelocity*zvelocity);
//    *tau    =  u_wall*(*xlength)/(C_S*C_S*Re)+0.5;
//    machNr  = u_wall/C_S;
//    printf("\nINFO: Calculated tau = %f \n", *tau);
//    printf("\nINFO: Wall speed = %f \n", u_wall);
//    printf("\nINFO: Mach number = %f \n\n", machNr);

//    /* valid settings check*/
//    if(u_wall >= C_S){
//    	ERROR("Wall speed is supersonic (aborting)! \n");
//    }
//
//    if(*tau<=0.5 || *tau>2){
//        ERROR("Tau is out of stability region (0.5,2.0) (aborting)! \n");
//    }
//
//    /*We allow user defined mach number tolerance for Ma << 1 (default = 0.1)
//      To change please look at LBDefinitions.h*/
//    if(machNr > machNrTol){
//        char buffer[80];
//        snprintf(buffer, 80, "Mach number is larger than %f (aborting)! \n",machNrTol);
//    	ERROR(buffer);
//    }

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

void initialiseFields(double *collideField, double *streamField, int *flagField,
		int *xlength, t_boundPara *boundPara, char *problem){


	// Loop variables
	int x,y,z;

    // Temporary variables for xlength^2
	int const xlen2		= xlength[0]+2;
	int const ylen2		= xlength[1]+2;
	int const zlen2		= xlength[2]+2;
	int offset1, offset2;

    int type;
    int start, end;

    /* TODO: (DL) there are 6x2 nested for loops now, to make it nicer
     * make a function and set the variables in the parameters.
     *
     * Problem is the different index computation in each type of wall -> use function pointer?
     */
    //Do the domain enclosing boundaries (ghost layers) first and set type accordingly.
    type = boundPara[XZ_BACK].type;
    start = boundPara[XZ_BACK].idxStartEnd[0];
    end = boundPara[XZ_BACK].idxStartEnd[1];
    for(z = start; z < zlen2+end; ++z){
		offset1 = z*xlen2*ylen2;
		offset2 = offset1 + (xlength[1]+1)*xlen2;
    	for(x = start; x < xlen2+end; ++x){
    		flagField[offset1 + x] = type; // y = 0
    	}
    }

    type = boundPara[XZ_FRONT].type;
    start = boundPara[XZ_FRONT].idxStartEnd[0];
    end = boundPara[XZ_FRONT].idxStartEnd[1];
    for(z = start; z < zlen2+end; ++z){
		offset1 = z*xlen2*ylen2;
		offset2 = offset1 + (xlength[1]+1)*xlen2;
    	for(x = start; x < xlen2+end; ++x){
    		flagField[offset2 + x] = type; // y = xlength[1]+1
    	}
    }

    type = boundPara[XY_LEFT].type;
    start = boundPara[XY_LEFT].idxStartEnd[0];
    end = boundPara[XY_LEFT].idxStartEnd[1];
    for(y = start; y < ylen2+end; ++y){
		offset1 = y*xlen2;
		offset2 = offset1 + (xlength[2]+1)*xlen2*ylen2;
    	for(x = start; x < xlen2+end; ++x){
    		flagField[offset1 + x] = type; // z = 0
    	}
    }

    type = boundPara[XY_RIGHT].type;
    start = boundPara[XY_RIGHT].idxStartEnd[0];
    end = boundPara[XY_RIGHT].idxStartEnd[1];
    for(y = start; y < ylen2+end; ++y){
		offset1 = y*xlen2;
		offset2 = offset1 + (xlength[2]+1)*xlen2*ylen2;
    	for(x = start; x < xlen2+end; ++x){
    		flagField[offset2 + x] = type; // z = xlength[2]+1
    	}
    }

    type = boundPara[YZ_BOTTOM].type;
    start = boundPara[YZ_BOTTOM].idxStartEnd[0];
    end = boundPara[YZ_BOTTOM].idxStartEnd[1];
    for(z = start; z < zlen2+end; ++z){
		offset1 = z*xlen2*ylen2;
		offset2 = offset1 + xlength[0] + 1;
    	for(y = start; y < ylen2+end; ++y){
			flagField[offset1 + y*xlen2] = type; // x = 0
    	}
    }

    type = boundPara[YZ_TOP].type;
    start = boundPara[YZ_TOP].idxStartEnd[0];
    end = boundPara[YZ_TOP].idxStartEnd[1];
    for(z = start; z < zlen2+end; ++z){
		offset1 = z*xlen2*ylen2;
		offset2 = offset1 + xlength[0] + 1;
    	for(y = start; y < ylen2+end; ++y){
    		flagField[offset2 + y*xlen2] = type; // x = xlength[0]+1
    	}
    }

    //NOTE: only domain (ghost layer were set previously)//
    for(z = 1; z <= xlength[2]; ++z){
		offset1 = z*xlen2*ylen2;
    	for (y = 1; y <= xlength[1]; ++y) {
			offset2 = offset1 + y*ylen2;
			for (x = 1; x <= xlength[0]; ++x) {
				int xyzoffset = offset2 + x;
				int type_domain = pgmMatrix[z][x];

				if(type_domain == 0){
					//This is actually not required as long as FLUID=0
					flagField[xyzoffset] = FLUID;
				}else if(type_domain == 1){
                    //TODO: Improve the checks for illegal geometry by saving these
					flagField[xyzoffset] = OBSTACLE;

				}else{
					ERROR("Description of scenario in pgm file should only consist of logical values. \n");
				}
			}
		}
    }

    print_flagfield_slice(flagField, xlength);

    /*Setting initial distributions, LATTICEWEIGHTS and INFLOW conditions */
    //f_i(x,0) = f^eq(1,0,0) = w_i

    // current cell index
    int idx;

    /* initialize collideField and streamField and inflow BC if any*/
    for ( z = 0; z < zlen2; ++z) {
		offset1 = z*xlen2*ylen2;
        for ( y = 0; y < ylen2; ++y) {
			offset2 = offset1 + y*xlen2;
            for ( x = 0; x < xlen2; ++x) {
				int xyzoffset =  offset2 + x;

				// Compute the base index
                idx = Q*xyzoffset;

                //Set initial condition
                if(flagField[xyzoffset] != INFLOW){
                    for (int i = 0; i < Q; ++i) {
                        collideField[idx+i] = LATTICEWEIGHTS[i];
                        streamField[idx+i]  = LATTICEWEIGHTS[i];
                    }
                }
                else{

                    /*Setting inflow condition once and for all*/
                    p_handleInflow(x, y, z, xlength, boundPara,
                        collideField,xyzoffset);
                }


                /* TODO: (DL) if we do not need p_handleInflow in boundary.c at all,
                 * we can also think to get it here.
                 */


            }
        }
    }

//	 for (x = 0; x < xlen2; x++) {
//         for (z = 0; z < zlen2; z++) {
//             printf("%d ",pgmMatrix[z][xlength[0]+1-x]);
//         }
//         printf("\n");
//     }
//
//	 printf("##############################################\n");
//
//	 //Compare with one layer in the flagfield:
//	 for(x = 0; x < xlen2; ++x){
//		 for(z = 0; z < zlen2; ++z){
//			 y = 3;
//			 printf("%d ",flagField[z*xlen2*ylen2 + y*xlen2 + x]);
//		 }
//		 printf("\n");
//	 }


    //Call at allocation
    //    pgmMatrix = imatrix(0, xlength[2]+1,0,xlength[0]+1);
    //    init_imatrix(pgmMatrix,0,xlength[2]+1,0,xlength[0]+1,1);
    //    print_matrix(xlength);
    free_imatrix(pgmMatrix,0,xlength[2]+1,0,xlength[0]+1);
}
