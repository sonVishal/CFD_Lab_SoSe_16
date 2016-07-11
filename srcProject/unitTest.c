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

void treatPostCollisionBoundary(const t_component * const c, t_procData const * const procData, int direction){
    t_iterPara iterPara1, iterPara2;
    int index1[5], index2[5];

    p_setCommIterationParameters(&iterPara1, procData, direction);
    p_setCommIterationParameters(&iterPara2, procData, direction+1);

    //opposite directions to indices
    p_assignIndices(direction+1,  index1);
    p_assignIndices(direction,    index2);

    int currentIndexField1 = -1, currentIndexField2=-1; //Initially assign to invalid

    static int shiftFixedValue[6] = {-1, 1, 1, -1, 1, -1}; //static to allocate/assign only once

    int newFixedValue1 = iterPara1.fixedValue + shiftFixedValue[direction];
    int newFixedValue2 = iterPara2.fixedValue + shiftFixedValue[direction+1];

    for(int m = 0; m < numComp; ++m){

        //k - corresponds to the 'outer' value when computing the offset
        for(int k = iterPara1.startOuter; k <= iterPara1.endOuter; ++k){
            //j - corresponds to the 'inner' value
            for(int j = iterPara1.startInner; j <= iterPara1.endInner; ++j){

                currentIndexField1  = Q*p_computeCellOffset(k, j, newFixedValue1, procData->xLength, direction);
                currentIndexField2  = Q*p_computeCellOffset(k, j, newFixedValue2, procData->xLength, direction+1);

                //out of bound checks
                assert(currentIndexField1 < Q*(procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2)
                        && currentIndexField1 >= 0);
                assert(currentIndexField2 < Q*(procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2)
                        && currentIndexField2 >= 0);

                double *currentCell1 = &c[m].streamField[currentIndexField1];
                double *currentCell2 = &c[m].streamField[currentIndexField2];

                for (int i = 0; i < 5; i++) {
                    currentCell1[index1[i]] = currentCell1[index1[i]] - (currentCell1[index1[i]] - c[m].feq[index1[i]])/(c[m].tau);
                    currentCell2[index2[i]] = currentCell2[index2[i]] - (currentCell2[index2[i]] - c[m].feq[index2[i]])/(c[m].tau);
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
        //Does also treat boundary for the respective counterpart (e.g. LEFT & RIGHT)
        treatPostCollisionBoundary(c, procData, LEFT);
        treatPostCollisionBoundary(c, procData, TOP);
        treatPostCollisionBoundary(c, procData, FRONT);
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
