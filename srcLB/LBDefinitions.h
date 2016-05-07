#include <math.h>

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

//TODO: (TKS) Consider #undef, but then we can't use them other places
//            in the code. Need to be careful with define.

#define w1 12.0/36.0
#define w2 2.0/36.0
#define w3 1.0/36.0

// Lattice weights
static const double LATTICEWEIGHTS[19] = {
    w3, w3, w2, w3, w3,
    w3, w2, w3, w2, w1,
    w2, w3, w2, w3, w3,
    w3, w2, w3, w3
};

// Speed of sound
//TODO: (TKS) For some reason the compiler complains that this is not a constant expression.
//static const double C_S = 1.0/sqrt(3.0);
// TODO: (VS) Remove once the above conflict is resolved
static const double C_S = 1.0;

#endif
