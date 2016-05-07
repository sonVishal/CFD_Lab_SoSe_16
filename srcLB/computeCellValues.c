#include "computeCellValues.h"

// TODO: (VS) Comments

// TODO: (VS) Find a clever way to avoid two functions as well as
    // conflicting data type error
void pDotProduct1(const double * const u, const double * const v,
    double * dotProd) {

    int i;
    *dotProd = 0.0;

    for (i = 0; i < 3; i++) {
        *dotProd += u[i]*v[i];
    }
}

void pDotProduct2(const double * const u, const int * const v,
    double * dotProd) {

    int i;
    *dotProd = 0.0;

    for (i = 0; i < 3; i++) {
        *dotProd += u[i]*v[i];
    }
}

void computeDensity(const double *const currentCell, double *density){
    const int Q = 19;
    int i;

    *density = 0.0;

    for (i = 0; i < Q; i++) {
        *density += currentCell[i];
    }
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
    int i;

    for (i = 0; i < Q; i++) {
        velocity[0] += currentCell[i]*LATTICEVELOCITIES[i][0];
        velocity[1] += currentCell[i]*LATTICEVELOCITIES[i][1];
        velocity[2] += currentCell[i]*LATTICEVELOCITIES[i][2];
    }

    velocity[0] /= *density;
    velocity[1] /= *density;
    velocity[2] /= *density;
}

void computeFeq(const double * const density, const double * const velocity, double *feq){
    int i;

    double dotProd1;
    double dotProd2;
    double C_S_2 = C_S*C_S;
    double tmp;

    *feq = 0.0;

    for (i = 0; i < Q; i++) {
        pDotProduct1(velocity, velocity, &dotProd1);
        pDotProduct2(velocity, LATTICEVELOCITIES[i], &dotProd2);
        tmp = dotProd2/C_S_2;
        *feq += LATTICEWEIGHTS[i]*(*density)*(1 + tmp + tmp*tmp/2 - dotProd1/2/C_S_2);
    }
}
