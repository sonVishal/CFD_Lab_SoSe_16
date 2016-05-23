#include "initLB.h"
#include "LBDefinitions.h"
#include "helper.h"

void p_readWall(char *argv[], t_boundPara *boundPara){
	int type;
	double x_velocity, y_velocity, z_velocity;
	double rhoRef, rhoIn;

    //TODO:(TKS) Only reads the first value named type. e.g. Only the first
    //           block of settings is red into the walls.
	READ_INT(*argv, type);

	READ_DOUBLE(*argv, x_velocity);
	READ_DOUBLE(*argv, y_velocity);
	READ_DOUBLE(*argv, z_velocity);

	READ_DOUBLE(*argv, rhoRef);
	READ_DOUBLE(*argv, rhoIn);

	boundPara->type = type;
	boundPara->wallVelocity[0] = x_velocity;
	boundPara->wallVelocity[1] = y_velocity,
	boundPara->wallVelocity[2] = z_velocity;
	boundPara->rhoRef = rhoRef;
	boundPara->rhoIn = rhoIn;

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

    double Re; //u_wall, machNr;
    int x_length, y_length, z_length; //Temporary for reading

    /* Read values from file given in argv */
    //Domain
    READ_INT(*argv, x_length);
    READ_INT(*argv, y_length);
    READ_INT(*argv, z_length);
    xlength[0] = x_length;
    xlength[1] = y_length;
    xlength[2] = z_length;

    //PGM-file that describes the scenario (located in /scenarios)
    READ_STRING(*argv, problem);

    for(int b=XY_LEFT; b <= XZ_BACK; b++){
    	p_readWall(argv, &boundPara[b]);
    }

    READ_DOUBLE(*argv, Re);
    READ_INT(*argv, *timesteps);
    READ_INT(*argv, *timestepsPerPlotting);

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

    /*Setting initial distributions*/
    //f_i(x,0) = f^eq(1,0,0) = w_i

    // current cell index
    int idx;

	// Loop variables
	int x,y,z;

    // Temporary variables for xlength^2
	int const xlen2		= xlength[0]+2;
	int const ylen2		= xlength[1]+2;
	int const zlen2		= xlength[2]+2;
	int offset1, offset2;

    /* initialize collideField and streamField */
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
                for (int i = 0; i < Q; ++i) {
                    collideField[idx+i] = LATTICEWEIGHTS[i];
                    streamField[idx+i]  = LATTICEWEIGHTS[i];
                }
            }
        }
    }

    int type1, type2;

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

    /* TODO: (DL) For the moment the pgm file is assumed to have the same
     * size as the values set in the parameter file (.dat).
     *
     * see issue: https://github.com/sonVishal/CFD_Lab_SoSe_16/issues/5
     *
     */
    int **pgmMatrix = read_pgm(problem);

    //NOTE: only domain (ghost layer were set previously)
    for(z = 1; z < zlen2; ++z){
		offset1 = z*xlen2*ylen2;
    	for (y = 1; y < ylen2; ++y) {
			offset2 = offset1 + y*ylen2;
			for (x = 1; x < xlen2; ++x) {
				int xyzoffset = offset2 + x;
				int type_domain = pgmMatrix[z][xlen2-1-x]; // values are mirrored
				flagField[xyzoffset] = type_domain;
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
