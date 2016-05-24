#include "initLB.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "boundary.h"
#include <stdio.h>

void p_readWall(char *argv[], t_boundPara *boundPara, const int skip){
	int type;
	double x_velocity, y_velocity, z_velocity;
	double rhoRef, rhoIn; /* TODO: (DL) rename rhoIn -> deltaRho, also in templates*/

    //TODO: Add a variable to the input saying how many variables of that should be skipped.
    //          * Remember to change both the #define part and the function definition.
    //          * Add if statement in all read functions.
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
}

//Checks wheter the sorroundings of the current cell is legal.
//Examples of geometries not allowed:
//  ##       #
//    ##     #
//            #
//            #

void valid_sorroundings(int *flagField, int currentCellIndex, int* xlength){

    int xlen2 = xlength[0] + 2;
    int ylen2 = xlength[1] + 2;


    //Check if YZ_TOPP edge has a FLUID neighbour.
    if(flagField[currentCellIndex + 1] == FLUID){

        //Check if XZ_BACK has a FLUID neighbour.
        if(flagField[currentCellIndex - xlen2] == FLUID){
            if(flagField[currentCellIndex + 1 - xlen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }

        //Check if XZ_FRONT has a FLUID neighbour.
        if(flagField[currentCellIndex + xlen2] == FLUID){
            if(flagField[currentCellIndex + 1 + xlen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }

        //Check if XY_LEFT has a FLUID neighbour.
        if(flagField[currentCellIndex - xlen2*ylen2] == FLUID){
            if(flagField[currentCellIndex + 1 - xlen2*ylen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }

        //Check if XY_RIGHT has a FLUID neighbour.
        if(flagField[currentCellIndex + xlen2*ylen2] == FLUID){
            if(flagField[currentCellIndex + 1 + xlen2*ylen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }

    }

    //Check if YZ_BOTTOM edge has a FLUID neighbour.
    if(flagField[currentCellIndex - 1] == FLUID){

        //Check if XZ_BACK has a FLUID neighbour.
        if(flagField[currentCellIndex - xlen2] == FLUID){
            if(flagField[currentCellIndex - 1 - xlen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }

        //Check if XZ_FRONT has a FLUID neighbour.
        if(flagField[currentCellIndex + xlen2] == FLUID){
            if(flagField[currentCellIndex - 1 + xlen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }

        //Check if XY_LEFT has a FLUID neighbour.
        if(flagField[currentCellIndex - xlen2*ylen2] == FLUID){
            if(flagField[currentCellIndex - 1 - xlen2*ylen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }

        //Check if XY_RIGHT has a FLUID neighbour.
        if(flagField[currentCellIndex + xlen2*ylen2] == FLUID){
            if(flagField[currentCellIndex - 1 + xlen2*ylen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }

    }

    //Check if XZ_FRONT has a FLUID neighbour
    if(flagField[currentCellIndex + xlen2] == FLUID){

        //Check if XY_RIGHT has a FLUID neighbour.
        if(flagField[currentCellIndex + xlen2*ylen2] == FLUID){
            if(flagField[currentCellIndex +xlen2 + xlen2*ylen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }

        //Check if XY_LEFT has a FLUID neighbour.
        if(flagField[currentCellIndex + xlen2*ylen2] == FLUID){
            if(flagField[currentCellIndex +xlen2 - xlen2*ylen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }
    }

    //Check if XZ_BACK has a FLUID neighbour
    if(flagField[currentCellIndex - xlen2] == FLUID){

        //Check if XY_RIGHT has a FLUID neighbour.
        if(flagField[currentCellIndex + xlen2*ylen2] == FLUID){
            if(flagField[currentCellIndex -xlen2 + xlen2*ylen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }

        //Check if XY_LEFT has a FLUID neighbour.
        if(flagField[currentCellIndex + xlen2*ylen2] == FLUID){
            if(flagField[currentCellIndex -xlen2 - xlen2*ylen2] == OBSTACLE)
                ERROR("Invalid geometry");
        }
    }


}

//Global variable makes it cleaner. A domain is described as in a pgm format
//even when no pgm file is provided. The two functions 'readParameters' and
//'initialiseFields' use it:
//'readParameters' - allocates and checks for validity of the geometry
//'initialiseFields' - uses the matrix to set the domain and deallocates the matrix
static int **pgmMatrix = NULL;

void init_pgmMatrix(const int * const xlength){

	pgmMatrix = imatrix(0, xlength[2]+2, 0, xlength[0]+2); //XZ-plane where Z is on the "x-axis"

	// set all values to 1 in the matrix
	init_imatrix(pgmMatrix, 0, xlength[2]+2, 0, xlength[0]+2, 1);

	// set all domain values to 0 (leaving the boundary to 1)
	init_imatrix(pgmMatrix, 1, xlength[2], 1, xlength[0], 0);
}

void read_customPgmMatrix(const int * const xlength, char *filename){

	char file_path[MAX_LINE_LENGTH+10]; //+10 for 'scenarios/'
	snprintf(file_path, MAX_LINE_LENGTH+10, "scenarios/%s", filename);
	int zsizePgm, xsizePgm;
	pgmMatrix = read_pgm(file_path, &zsizePgm, &xsizePgm);

	if(zsizePgm != xlength[2] && xsizePgm != xlength[0]){
		free_imatrix(pgmMatrix,0,zsizePgm+2,0,xsizePgm+2);
		ERROR("TODO: msg");
	}

	/* TODO: (DL) DO CHECK FOR VALID GEMOETRY HERE
	 * Also remember that we need to free the matrix, so it's better to
	 * return an error code..
	 */
}

/* TODO: (DL) delete when finializing code, but leave it as long debugging there
* is possible debugging going on :) */
void print_matrix(const int * const xlength){
	 for (int x = 0; x < xlength[0]+2; x++) {
       for (int z = 0; z < xlength[2]+2; z++) {
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
			//TODO: (DL) make error msg with values..
			ERROR("The step cannot be bigger than the domain.");
		}

		init_imatrix(pgmMatrix, 1, z_direction, 1, x_direction, OBSTACLE);
		break;
	case ARBITRARY:
		READ_STRING(*argv, problem, 0);
		read_customPgmMatrix(xlength, problem);
		break;
	default:
		ERROR("SETTING-MODE is not known!");
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
    if(MODE < 0 || MODE > NUM_MODES ){
    	char msg[40];
    	snprintf(msg, 40, "SETTING-MODE=%i is not known! \n", MODE);
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

void initialiseFields(double *collideField, double *streamField, int *flagField,
		int *xlength, t_boundPara *boundPara, char *problem){


	// Loop variables
	int x,y,z;

    // Temporary variables for xlength^2
	int const xlen2		= xlength[0]+2;
	int const ylen2		= xlength[1]+2;
	int const zlen2		= xlength[2]+2;
	int offset1, offset2;

    int type1, type2;

    /* TODO: (DL) How to deal with overlapping cells? There are cells at edges and corners that
     * are included twice. The final value these cells is determined by the last loop that
     * contains these cells.
     */

    //Do the domain enclosing boundaries (ghost layers) first and set type accordingly.
    type1 = boundPara[XZ_FRONT].type;
    type2 = boundPara[XZ_BACK].type;
    for(z = 0; z < zlen2; ++z){
		offset1 = z*xlen2*ylen2;
		offset2 = offset1 + (xlength[1]+1)*xlen2;
    	for(x = 0; x < xlen2; ++x){
    		flagField[offset1 + x] = type1; // y = 0
    		flagField[offset2 + x] = type2; // y = xlength[1]+1
    	}
    }

    type1 = boundPara[XY_LEFT].type;
    type2 = boundPara[XY_RIGHT].type;
    for(y = 0; y < ylen2; ++y){
		offset1 = y*xlen2;
		offset2 = offset1 + (xlength[2]+1)*xlen2*ylen2;
    	for(x = 0; x < xlen2; ++x){
    		flagField[offset1 + x] = type1; // z = 0
    		flagField[offset2 + x] = type2; // z = xlength[2]+1
    	}
    }

    type1 = boundPara[YZ_BOTTOM].type;
    type2 = boundPara[YZ_TOP].type;
    for(z = 0; z < zlen2; ++z){
		offset1 = z*xlen2*ylen2;
		offset2 = offset1 + xlength[0] + 1;
    	for(y = 0; y < ylen2; ++y){
			flagField[offset1 + y*xlen2] = type1; // x = 0
    		flagField[offset2 + y*xlen2] = type2; // x = xlength[0]+1
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

    /*Check for illegal geometries*/
    // TODO: (TKS) Could do this in the above loop be checking pgmMatrix

    for(z = 1; z <= xlength[2]; ++z){
		offset1 = z*xlen2*ylen2;
    	for (y = 1; y <= xlength[1]; ++y) {
			offset2 = offset1 + y*ylen2;
			for (x = 1; x <= xlength[0]; ++x) {
				int xyzoffset = offset2 + x;
				if(flagField[xyzoffset] == OBSTACLE){
                    valid_sorroundings(flagField, xyzoffset, xlength);
				}
			}
		}
    }

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
				/*TODO: (DL) maybe it's required to set certain boundaries/obstacles
				 * differently... in the cavity we set also the boundaries, so for now
				 * there is no check
				 */
				// Compute the base index
                idx = Q*xyzoffset;

                //Set initial condition
                for (int i = 0; i < Q; ++i) {
                    collideField[idx+i] = LATTICEWEIGHTS[i];
                    streamField[idx+i]  = LATTICEWEIGHTS[i];
                }

                /* TODO: (DL) not sure if that was supposed to be like this?
                 * But can the below if-statement decide which one of the two
                 * loops to execute. At the moment it seems it first gets set to
                 * LATTICEWEIGHTS and then to inflow condition (if it's an INFLOW).
                 */

                /* TODO: (DL) if we do not need p_handleInflow in boundary.c at all,
                 * we can also think to get it here.
                 */

                /*Setting inflow condition once and for all*/
                if(flagField[idx] == INFLOW){
                    for (int i = 0; i < Q; ++i) {
                        p_handleInflow( x,  y,  z,  xlength, boundPara, 
                                        collideField, idx);
                    }
                }

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

    free_imatrix(pgmMatrix,0,zlen2,0,xlen2);
}
