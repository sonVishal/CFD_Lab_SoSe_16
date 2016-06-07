#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

// Number of lattice directions
static const int Q = 19;

// Lattice velocity vectors
static const int LATTICEVELOCITIES[19][3] = {
    {0,-1,-1}, {-1,0,-1}, {0,0,-1}, {1,0,-1}, {0,1,-1},
    {-1,-1,0}, {0,-1,0},  {1,-1,0}, {-1,0,0}, {0,0,0},
    {1,0,0},   {-1,1,0},  {0,1,0},  {1,1,0},  {0,-1,1},
    {-1,0,1},  {0,0,1},   {1,0,1},  {0,1,1}
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
enum {
    LEFT,
    RIGHT,
    TOP,
    BOTTOM,
    FRONT,
    BACK
};

// Enum for cell type
enum {
    FLUID,
    NO_SLIP,
    MOVING_WALL,
    PARALLEL_BOUNDARY
};

// Local struct for each process
typedef struct {
    int rank;
    int numRanks;
    int xLength[3];
    int neighbours[6]; //same order as enum
    double wallVelocity[3];
} t_procData;

// Function for reverse index search for procs
static inline void p_indexToPos(const int *const len, const int myIndex, int *myPos) {
    int i,j,k;
    for (k = 0; k < len[2]; k++) {
        for (j = 0; j < len[1]; j++) {
            for (i = 0; i < len[0]; i++) {
                if (i+j*len[0]+k*len[0]*len[1] == myIndex) {
                    myPos[0] = i;
                    myPos[1] = j;
                    myPos[2] = k;
                    return;
                }
            }
        }
    }
}
#endif
