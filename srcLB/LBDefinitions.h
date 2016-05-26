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

//Note: better not change the order, because at the moment it is defined in the parameter
//files using the integer value!
enum E_CELL_TYPE{
	FLUID,
	NO_SLIP,
	MOVING,
	FREE_SLIP,
	INFLOW,
	OUTFLOW,
	PRESSURE_IN,
	OBSTACLE
};

//This order of boundaries has to correspond to the order in the parameter file (.dat file)!!
enum E_BOUNDARIES{
	XY_LEFT,
	XY_RIGHT,
	YZ_BOTTOM,
	YZ_TOP,
	XZ_FRONT,
	XZ_BACK,
	NUM_WALLS //would be currently 6 (to have this dynamically, can be used in for loops)
};

typedef struct {
	int type;
	int idxStartEnd[2];
    // Three first indecies hold the velocities in each direction and the fourth the magnitude
	double wallVelocity[3];
    double rhoRef;
	double rhoIn;
} t_boundPara;

typedef struct {
    short type;
    short position;
} t_flagField;

enum SETTING_MODE{
	SHEAR_FLOW,
	STEP_FLOW,
	ARBITRARY, 	//a .pgm file has to be provided in this mode (e.g. example tilted plate)
    CAVITY,
	NUM_MODES //used for inbound checks - if inserting a new mode, that just before this one!
};

// Typedef for the function pointer used in treating boundaries
typedef void (*t_boundaryFcnPtr)(double*, t_flagField const * const, int const * const,
                int const * const, const t_boundPara * const,
                int const * const);


// Functions for calculating indecies
inline static void p_computeIndex(const int * const point,
	const int * const xlength, int *index) {
	(*index) = (xlength[0]+2)*((xlength[1]+2)*point[2]+point[1]) + point[0];
}

inline static void p_computeIndexXYZ(const int x, const int y, const int z,
	const int * const xlength, int *index) {
	(*index) = (xlength[0]+2)*((xlength[1]+2)*z+y) + x;
}

inline static void p_computeIndexQ(const int * const point,
	const int * const xlength, int *index) {
	(*index) = Q*((xlength[0]+2)*((xlength[1]+2)*point[2]+point[1]) + point[0]);
}
#endif
