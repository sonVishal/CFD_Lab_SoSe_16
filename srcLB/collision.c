#include "collision.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
    int i;
    for (i = 0; i < 19; i++) {
        currentCell[i] = currentCell[i] - (currentCell[i]-feq[i])/(*tau);
    }
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){
    /* TODO */
}
