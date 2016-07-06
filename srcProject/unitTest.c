#include "unitTest.h"

double **massBefore;
double **massAfter;
double momentumBefore[3];
double momentumAfter[3];

void initializeUnitTest(const int totalSize){
    massBefore = (double **) malloc(numComp * sizeof(double*));
    massAfter = (double **) malloc(numComp * sizeof(double*));

    for (int i = 0; i < numComp; i++) {
        massBefore[i] = (double *)  malloc(totalSize * sizeof( double ));
        massAfter[i] = (double *)  malloc(totalSize * sizeof( double ));
    }
}

void storeMassVector(const t_component * const c, const int * const xlength, double **massVector) {
    // Begin
    int idx;
    for (int i = 0; i < numComp; i++) {
        for (int z = 1; z <= xlength[2]; z++) {
            for (int y = 1; y <= xlength[1]; y++) {
                for (int x = 1; x <= xlength[0]; x++) {
                    idx = p_computeCellOffsetXYZ(x, y, z, xlength);
                    massVector[i][idx] = c[i].rho[idx];
                }
            }
        }
    }
}

void computePostCollisionDistributions(const t_component * const c, t_procData const * const procData){

    for(int k = 0; k < numComp; ++k){

        for (int z = 1; z <= procData->xLength[2] ; z++) {
        	for (int y = 1; y <= procData->xLength[1]; y++) {
        		for (int x = 1; x <= procData->xLength[0]; x++) {

                    int idx = p_computeCellOffsetXYZ_Q(x, y, z, procData->xLength);
                    double *currentCell = &c[k].streamField[idx];

                    for(int i=0; i<Q ; ++i){
            		    currentCell[i] = currentCell[i] - (currentCell[i] - c[k].feq[i])/(c[k].tau);
                        assert(currentCell[i] >= 0);
                    }

            	}
            }
        }
    }

}

void doStreaming(const t_component * const c, t_procData const * const procData){
    int nextCellIndex, cellIdx;

    for(int k=0; k < numComp; ++k){

        for (int z = 1; z <= procData->xLength[2] ; z++) {
        	for (int y = 1; y <= procData->xLength[1]; y++) {
        		for (int x = 1; x <= procData->xLength[0]; x++) {

                    cellIdx = p_computeCellOffsetXYZ_Q(x, y, z, procData->xLength);

                    for(int j=0; j<Q ; ++j){
                        nextCellIndex = p_computeCellOffsetXYZ_Q(x+LATTICEVELOCITIES[j][0],
                            y+LATTICEVELOCITIES[j][1], z+LATTICEVELOCITIES[j][2], procData->xLength);

                        c[k].collideField[cellIdx+Q-j-1] = c[k].streamField[nextCellIndex+Q-j-1];
                    }

            	}
            }
        }
    }
}

void steamCollideUnitTest(const t_component * const c, t_procData const * const procData){

    beforeCollision(c, procData);
    computePostCollisionDistributions(c, procData);
    afterCollision(c, procData);

    doStreaming(c, procData);
}

void checkMassVector(const int * const xlength, const int rank) {
    // Begin
    int idx;
    for (int i = 0; i < numComp; i++) {
        for (int z = 1; z <= xlength[2]; z++) {
            for (int y = 1; y <= xlength[1]; y++) {
                for (int x = 1; x <= xlength[0]; x++) {
                    idx = p_computeCellOffsetXYZ(x, y, z, xlength);
                    if (fabs(massAfter[i][idx] - massBefore[i][idx]) > _TOL_) {
                        printf("Proc: %d, Mass not conserved in cell (%d,%d,%d)\n",rank,x,y,z);
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

void beforeCollision(const t_component * const c, t_procData const * const procData){
    storeMassVector(c, procData->xLength, massBefore);
    computeGlobalMomentum(c, procData->xLength, momentumBefore);
}

void afterCollision(const t_component * const c, t_procData const * const procData){
    storeMassVector(c, procData->xLength, massAfter);
    computeGlobalMomentum(c, procData->xLength, momentumAfter);

    checkMassVector(procData->xLength, procData->rank);
    if (procData->rank == 0) {
        checkMomentum();
    }
}

void computeGlobalMomentum(const t_component *const c, const int *const xlength, double *compMomentum){
    //Begin
    int idx;
    double tempMomentum[3] = {0.0};

    for (int i = 0; i < numComp; i++) {
        for (int z = 1; z <= xlength[2]; z++) {
            for (int y = 1; y <= xlength[1]; y++) {
                for (int x = 1; x <= xlength[0]; x++) {
                    idx = p_computeCellOffsetXYZ_Q(x, y, z, xlength);
                    double cellMomentum[3];

                    //TODO: (check if streamField is correct)...
                    computeCellMomentum(&c[i].streamField[idx], cellMomentum);
                    tempMomentum[0] += c[i].m*cellMomentum[0];
                    tempMomentum[1] += c[i].m*cellMomentum[1];
                    tempMomentum[2] += c[i].m*cellMomentum[2];
                }
            }
        }
    }
    MPI_Reduce(tempMomentum, compMomentum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

void checkMomentum() {
    for (int i = 0; i < 3; i++) {
        if (fabs(momentumAfter[i]-momentumBefore[i]) > _TOL_) {
            printf("Global momentum not conserved\n");
        }
    }
}

void freeUnitTest(){
    for (int i = 0; i < numComp; i++) {
        free(massBefore[i]);
        free(massAfter[i]);
    }
    free(massBefore);
    free(massAfter);
}
