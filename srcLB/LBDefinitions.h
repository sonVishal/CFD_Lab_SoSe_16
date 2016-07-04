#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_
#include <math.h>

// Number of lattice directions
static const int Q = 19;

// reference density
static const double rhoRef = 0.693;
static const int rho0   = 1;

#define NUMCOMP 1

// Lattice velocity vectors
static const int LATTICEVELOCITIES[19][3] = {
    {0,-1,-1}, {-1,0,-1}, {0,0,-1}, {1,0,-1}, {0,1,-1},
    {-1,-1,0}, {0,-1,0},  {1,-1,0}, {-1,0,0}, {0,0,0},
    {1,0,0},   {-1,1,0},  {0,1,0},  {1,1,0},  {0,-1,1},
    {-1,0,1},  {0,0,1},   {1,0,1},  {0,1,1}
};

// Temporary variables for Lattice weights
#define w1 1.0/3.0  /* if ||c_i|| = 0 */
#define w2 1.0/18.0   /* if ||c_i|| = 1 */
#define w3 1.0/36.0   /* if ||c_i|| = sqrt(2) */

// Lattice weights
static const double LATTICEWEIGHTS[19] = {
    w3, w3, w2, w3, w3,
    w3, w2, w3, w2, w1,
    w2, w3, w2, w3, w3,
    w3, w2, w3, w3
};


typedef struct {
    double* streamField;
    double* collideField;
    double* feq;
    double* rho;
    double** velocity;
    double** force;
    double  tau;
    double  m;
    double  d0;
    int psiFctCode;
} t_component;

//TODO: For now assuming that mass = 1. The input density if mass changes.
static inline double psi0(double numberDensity){
    // const double E = 2.71828;
    const double rho0 = 1.0;
    // return rho0*(1 - pow(E, -numberDensity/rho0));
    return rho0*(1-exp(-numberDensity/rho0));
}
static inline double psi1(double numberDensity){ return numberDensity;}
static inline double psi2(double numberDensity){ return 4.0*exp(-rho0/numberDensity);}

typedef double (*fctPtrPsi)(double);
static const fctPtrPsi psiFctPointer[3] = {psi0, psi1, psi2};

static inline int p_computeCellOffsetXYZ_Q(const int x, const int y, const int z, const int xlength) {
    return Q*(x+(y+z*(xlength+2))*(xlength+2));
}

static inline int p_computeCellOffsetXYZ(const int x, const int y, const int z, const int xlength) {
    return x+(y+z*(xlength+2))*(xlength+2);
}

enum {
    FLUID,
    MOVING_WALL,
    NO_SLIP,
    PERIODIC
};

enum {
    LEFT,
    RIGHT,
    TOP,
    BOTTOM,
    FRONT,
    BACK
};

// Speed of sound
#define C_S ((double)0.57735026918962576) // double cast to guarantee type safety

// Tolerances
static const double densityTol  = 0.03;
static const double machNrTol   = 0.1;

#endif
