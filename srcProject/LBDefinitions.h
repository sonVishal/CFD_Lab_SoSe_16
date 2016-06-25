#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_
#include "helper.h"

// Number of lattice directions
static const int Q = 19;

// Number of components is a global variable
int numComp;

// Lattice velocity vectors
static const int LATTICEVELOCITIES[19][3] = {
    {0,-1,-1}, {-1,0,-1}, {0,0,-1}, {1,0,-1}, {0,1,-1}, //[0-4]
    {-1,-1,0}, {0,-1,0},  {1,-1,0}, {-1,0,0}, {0,0,0},  //[5-9]
    {1,0,0},   {-1,1,0},  {0,1,0},  {1,1,0},  {0,-1,1}, //[10-14]
    {-1,0,1},  {0,0,1},   {1,0,1},  {0,1,1}             //[15-18]
};


// Temporary variables for Lattice weights
#define w1 12.0/36.0  /* if ||c_i|| = 0 */
#define w2 2.0/36.0   /* if ||c_i|| = 1 */
#define w3 1.0/36.0   /* if ||c_i|| = sqrt(2) */

// Lattice weights
static const double LATTICEWEIGHTS[19] = {
    w3, w3, w2, w3, w3,
    w3, w2, w3, w2, w1,
    w2, w3, w2, w3, w3,
    w3, w2, w3, w3
};

#undef w1
#undef w2
#undef w3

// Speed of sound
#define C_S ((double)0.57735026918962576) // double cast to guarantee type safety

// Tolerances
static const double densityTol  = 0.03;
static const double machNrTol   = 0.1;

// Enum for send and read buffer directions
enum WALLS {
    LEFT,
    RIGHT,
    TOP,
    BOTTOM,
    FRONT,
    BACK
};

// Enum for cell type
enum CELLS {
    FLUID,
    NO_SLIP,
    MOVING_WALL,
    PARALLEL_BOUNDARY,
    PERIODIC_BOUNDARY
};

//TODO: (DL) not sure if a function pointer with inline works... we can test it out.
static inline double psi1(double numberDensity) {
    return (1-exp(-numberDensity));
}

static inline double psi2(double numberDensity) {
    return numberDensity;
}

typedef double (*fctPtrPsi)(double);
fctPtrPsi selectPsiFunction(int code){
    switch(code){
        case 1:
            return &psi1;
        case 2:
            return &psi2;
        default:
            ERROR("Function code for psi is invalid.");
    }
    return NULL;
}


//TODO: (DL) There can be "better" ways to number edges. At the moment the numbering is benefical for computing indices.
// Numbering of edges of the domain:
// //
//
//            ^  (Z)
//            |
//            |
//            |
//            |             (11)
//            +-----------------------------+
//        (8)/|                            /|
//          / |        (9)            (10)/ |
//         +-----------------------------+  |
//         |  |                          |  |
//         |  | (4)                      |  |(7)
//         |  |                          |  |
//      (5)|  |                       (6)|  |                    (Y)
//         |  +--------------------------|--+--------------------->
//         | /           (3)             | /
//         |/(0)                         |/ (2)
//         +-----------------------------+
//        /          (1)
//       / (X)
//      v

// Local struct for each process
typedef struct {
    int rank;
    int numRanks;
    int xLength[3];
    int neighbours[6];         //same order as enum //TODO: (DL) possibly rename neighbours, e.g. "parallelNeighbours"
    int periodicNeighbours[6]; //only value if MPI_PROC_NULL in neighbours
    int periodicEdgeNeighbours[12]; //shared edges - diagonal opposite neighbours
    int bufferSize[3]; //buffer sizes for all directions
    double wallVelocity[3];
} t_procData;

typedef struct {
    double* streamField;
    double* collideField;
    double  tau;
    double  m;
    double (*psi)(double);
} t_component;

/* Function for reverse index search for procs.
 * Is required in parallel.c and visualLB.c
 */
static inline void p_rankToPos(const int *const len, const int rank, int *pos) {
    int i,j,k;
    for (k = 0; k < len[2]; k++) {
        for (j = 0; j < len[1]; j++) {
            for (i = 0; i < len[0]; i++) {
                if (i+j*len[0]+k*len[0]*len[1] == rank) {
                    pos[0] = i;
                    pos[1] = j;
                    pos[2] = k;
                    return;
                }
            }
        }
    }
}

static inline void p_setIterationParameters(int *endOuter, int *endInner, int *fixedValue, const t_procData * const procData, const int wallIdx){

	switch(wallIdx/2){
	//---------------------------------------------
	//outer = Z, inner = X, Y fixed
	case 0:
		*endOuter = procData->xLength[2]+1;
		*endInner = procData->xLength[0]+1;
		*fixedValue = wallIdx == LEFT ? 0 : procData->xLength[1]+1;
		break;

	//---------------------------------------------
	//outer = Y, inner = X, Z fixed
	case 1:
		*endOuter = procData->xLength[1]+1;
		*endInner = procData->xLength[0]+1;
		*fixedValue = wallIdx == BOTTOM ? 0 : procData->xLength[2]+1;
		break;

	//---------------------------------------------
	//outer = Z, inner = Y, X fixed
	case 2:
		*endOuter = procData->xLength[2]+1;
		*endInner = procData->xLength[1]+1;
		*fixedValue = wallIdx == BACK ? 0 : procData->xLength[0]+1;
		break;

    default:
	// 	ERROR("Invalid wallIdx occurred. This should never happen!");
        assert(wallIdx <= 6 && wallIdx >= 0);
    }
}

/* Helper function that computes the offset of the current cell. 'Inner' corresponds to the first value of (x,y,z)
 * that is not fixed; 'outer' to the second value. By this ordering a better cache efficiency is obtained.
 * WallIdx has to be a valid index from the enum WALLS.
 */
static inline int p_computeCellOffset(const int outer, const int inner, const int fixedValue, int const * const xlength, const int wallIdx){

	//Computes index: z * (xlen*ylen) + y * (xlen) + x

	//wallIdx has valid integer values from 0 to 5
	switch (wallIdx/2) { //integer division to get the type of face (see enum in LBDefinitions.h)
		case 0: // LEFT, RIGHT -> Y fixed
			//outer = Z, inner = X
			return (xlength[0]+2) * (outer * (xlength[1]+2) + fixedValue) + inner;

		case 1: // TOP, BOTTOM -> Z fixed
			//outer = Y, inner = X
			return (xlength[0]+2) * (fixedValue * (xlength[1]+2) + outer) + inner;

		case 2: // FRONT, BACK -> X fixed
			//outer = Z, inner = Y
			return (xlength[0]+2) * (outer * (xlength[1]+2) + inner) + fixedValue;

		default:
			// ERROR("Invalid wall index occured. This should not happen !!!");
            assert(wallIdx <= 6 && wallIdx >= 0);
			return -1;
	}
}

static inline int p_computeCellOffsetXYZ(const int x, const int y, const int z, const int * const xlength) {
    return Q*(x+(y+z*(xlength[1]+2))*(xlength[0]+2));
}

#endif
