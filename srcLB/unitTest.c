#include "unitTest.h"
#include <stdlib.h>

double **massBefore;
double **massAfter;
double momentumBefore[3];
double momentumAfter[3];


void initializeUnitTest(const int totalSize){
    massBefore = (double **) malloc(1 * sizeof(double*));
    massAfter = (double **) malloc(1 * sizeof(double*));

    for (int i = 0; i < 1; i++) {
        massBefore[i] = (double *)  malloc(totalSize * sizeof( double ));
        massAfter[i] = (double *)  malloc(totalSize * sizeof( double ));
    }
}

void storeMassVector(const t_component * const c, const int xlength, double **massVector) {
    // Begin
    int idx;
    for (int i = 0; i < 1; i++) {
        for (int z = 1; z <= xlength; z++) {
            for (int y = 1; y <= xlength; y++) {
                for (int x = 1; x <= xlength; x++) {
                    idx = p_computeCellOffsetXYZ(x, y, z, xlength);
                    computeNumDensity(&c[i].collideField[Q*idx], &massVector[i][idx]);
                }
            }
        }
    }
}

void checkMassVector(const int xlength) {
    // Begin
    int idx;
    for (int i = 0; i < 1; i++) {
        for (int z = 1; z <= xlength; z++) {
            for (int y = 1; y <= xlength; y++) {
                for (int x = 1; x <= xlength; x++) {
                    idx = p_computeCellOffsetXYZ(x, y, z, xlength);
                    if (fabs(massAfter[i][idx] - massBefore[i][idx]) > _TOL_) {
                        printf("Proc: Mass not conserved in cell (%d,%d,%d)\n", x,y,z);
                    }
                }
            }
        }
    }
}

void computeCellMomentum(const double * const currentCell, double *momentum) {
    // Compute the momentum as sum(f_i*c_i)
    momentum[0] = -currentCell[1]+currentCell[3];
    momentum[1] = -currentCell[0];
    momentum[2] = -(currentCell[0]+currentCell[1]+currentCell[2]+currentCell[3]);

    // 4 to 7
    momentum[0] += -currentCell[5]+currentCell[7];
    momentum[1] += currentCell[4]-(currentCell[5]+currentCell[6]+currentCell[7]);
    momentum[2] += -currentCell[4];

    // 8 to 11
    momentum[0] += -currentCell[8]+currentCell[10]-currentCell[11];
    momentum[1] += currentCell[11];

    // 12 to 15
    momentum[0] += currentCell[13]-currentCell[15];
    momentum[1] += currentCell[12]+currentCell[13]-currentCell[14];
    momentum[2] += currentCell[14]+currentCell[15];

    // 16 to 18
    momentum[0] += currentCell[17];
    momentum[1] += currentCell[18];
    momentum[2] += currentCell[16]+currentCell[17]+currentCell[18];
}

void beforeCollision(const t_component * const c, int xLength){
    storeMassVector(c, xLength, massBefore);
    computeGlobalMomentum(c, xLength, momentumBefore);
}

void afterCollision(const t_component * const c,const int xLength){
    storeMassVector(c, xLength, massAfter);
    computeGlobalMomentum(c, xLength, momentumAfter);
    checkMassVector(xLength);
    checkMomentum();
}

void computeGlobalMomentum(const t_component *const c, const int xlength, double *compMomentum){
    // Begin
    int idx;
    double tempMomentum[3] = {0.0};

    int i = 0;

    for (int z = 1; z <= xlength; z++) {
        for (int y = 1; y <= xlength; y++) {
            for (int x = 1; x <= xlength; x++) {
                idx = p_computeCellOffsetXYZ_Q(x, y, z, xlength);
                double cellMomentum[3];
                computeCellMomentum(&c[i].collideField[idx], cellMomentum);
                tempMomentum[0] += c[i].m*cellMomentum[0];
                tempMomentum[1] += c[i].m*cellMomentum[1];
                tempMomentum[2] += c[i].m*cellMomentum[2];
            }
        }
    }
}

void checkMomentum() {
    for (int i = 0; i < 3; i++) {
        if (fabs(momentumAfter[i]-momentumBefore[i]) > _TOL_) {
            printf("Global momentum not conserved\n");
        }
    }
}

void freeUnitTest(){
    for (int i = 0; i < 1; i++) {
        free(massBefore[i]);
        free(massAfter[i]);
    }
    free(massBefore);
    free(massAfter);
}
