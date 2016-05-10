#include "computeCellValues.h"

// TODO: (VS) Find a clever way to avoid two functions as well as
    // conflicting data type error

// Function to compute dot product of 3D vectors u and v
void pDotProduct1(const double * const u, const double * const v,
    double * dotProd) {

    int i;
    *dotProd = 0.0;

    for (i = 0; i < 3; i++) {
        *dotProd += u[i]*v[i];
    }
}

// Function to compute dot product of 3D vectors u and lattice velocities c[i]
void pDotProduct2(const double * const u, const int * const v,
    double * dotProd) {

    int i;
    *dotProd = 0.0;

    for (i = 0; i < 3; i++) {
        *dotProd += u[i]*v[i];
    }
}

/** computes the density from the particle distribution functions stored at
 *  currentCell. currentCell thus denotes the address of the first particle
 *  distribution function of the respective cell.
 *  The result is stored in density.
 */
void computeDensity(const double *const currentCell, double *density){
    int i;

    // Initialize to 0.0
    *density = 0.0;

    // Compute the cell density
    for (i = 0; i < Q; i++) {
        *density += currentCell[i];
    }
}

/** computes the velocity within currentCell and stores the result in velocity */
void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
    int i;

    // Compute the cell momentum
    for (i = 0; i < Q; i++) {
        velocity[0] += currentCell[i]*LATTICEVELOCITIES[i][0];
        velocity[1] += currentCell[i]*LATTICEVELOCITIES[i][1];
        velocity[2] += currentCell[i]*LATTICEVELOCITIES[i][2];
    }

    // Divide the momentum by density to get velocities
    velocity[0] /= *density;
    velocity[1] /= *density;
    velocity[2] /= *density;
}

/** computes the equilibrium distributions for all particle distribution
 *  functions of one cell from density and velocity and stores the results in feq.
 */
void computeFeq(const double * const density, const double * const velocity, double *feq){
    int i;

    // Temporary variables
    double dotProd1;
    double dotProd2;
    double C_S_2 = C_S*C_S;
    double tmp;

    // Compute the equilibrium distribution for each component
    for (i = 0; i < Q; i++) {
        pDotProduct1(velocity, velocity, &dotProd1);
        pDotProduct2(velocity, LATTICEVELOCITIES[i], &dotProd2);
        tmp = dotProd2/C_S_2;
        feq[i] = LATTICEWEIGHTS[i]*(*density)*(1 + tmp + tmp*tmp/2 - dotProd1/2/C_S_2);
    }
}
