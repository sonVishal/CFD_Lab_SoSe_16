
#include "LBDefinitions.h"
#include "helperMPI.h"
#include "boundary.h"
#include "helper.h"
#include <mpi/mpi.h>

//TODO: (TKS) 
//Reduce bufferSize[6]--> bufferSize[3]?
//          * Can do two iterations of the loop at the same time.
//          * Or complete loop unroll
//          * REMOVE it if not used
//
//At what level to iterate plane?
//          * inside each function
//          * make a double for in communicate and use inner/outer formulation.
//                  * Extract/Inject need only get indecies of cillideField + buffer.
//


void communicate(double** sendBuffer, double**readBuffer, double* collideField, const t_procData *procData){
    //Run extract, swap, inject for all sides and cells.
    //
    int index[5];

    t_iterPara iterPara;


    for (int direction = LEFT; direction <= BACK; direction+=2) {
        p_setCommIterationParameters(&iterPara, procData, direction);
        p_assignIndices(&direction, index);

        extract(sendBuffer, collideField, &iterPara, procData, direction, index);
        extract(sendBuffer, collideField, &iterPara, procData, direction+1, index);
        swap(sendBuffer, readBuffer, procData, &direction);
        inject(readBuffer, collideField, &iterPara, procData, direction, index);
        inject(readBuffer, collideField, &iterPara, procData, direction+1, index);
    }

}

//Copy distributions needed to sendbuffer.
void extract( double** sendBuffer, double* collideField, const t_iterPara *iterPara, const t_procData *procData,
              int direction, int* index){

    int currentIndexField;
    int currentIndexBuff;

    //For error checking. may remove later.
    int fieldSize = Q*(procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2);
    //printf("fieldSize = %d\n", fieldSize);

    //k - corresponds to the 'outer' value when computing the offset
    for(int k = iterPara->startOuter; k <= iterPara->endOuter; ++k){
        //j - corresponds to the 'inner' value
        for(int j = iterPara->startInner; j <= iterPara->endInner; ++j){

            currentIndexField  = Q*p_computeCellOffset(k, j, iterPara->fixedValue, procData->xLength, direction);
            currentIndexBuff  =  5*p_computeBuffCellOffset(k, j, procData->bufferLength, direction);

            assert(currentIndexBuff < procData->bufferSize[direction/2] && currentIndexBuff >=0);
            //printf("currentIndexField = %d\n", currentIndexField);
            assert(currentIndexField < fieldSize  && currentIndexField <= 0);

            for (int i = 0; i < 5; i++) {
                //TODO: (TKS) Not correct indexing yet
                sendBuffer[direction][currentIndexBuff + i] = collideField[currentIndexField+index[i]];
            }
        }
    }
}

//Send distributions and wait to receive.
void swap(double** sendBuffer, double**readBuffer, const t_procData *procData, int *direction){

}

//Copy read buffer to ghost layer
void inject(double** readBuffer, double* collideField, const t_iterPara *iterPara, const t_procData *procData, 
            int direction, int *index){

    //k - corresponds to the 'outer' value when computing the offset
    for(int k = iterPara->startOuter; k <= iterPara->endOuter; ++k){
        //j - corresponds to the 'inner' value
        for(int j = iterPara->startInner; j <= iterPara->endInner; ++j){

          //int currentCellIndex = Q*p_computeCellOffset(k, j, iterPara.fixedValue, procData->xLength, direction);

        }
    }
}

//Function to assign iteration parameters for communication.
void p_setCommIterationParameters(t_iterPara *iterPara, const t_procData *procData, const int direction){
    //TODO: (TKS) Confirm that indecies are correct.
	switch(direction){

	//---------------------------------------------
	//outer = Z, inner = X, Y fixed
    //only iterate over inner domain of plane (FLUID cells)
	case LEFT:
        iterPara->startOuter = 1;
		iterPara->endOuter   = procData->xLength[2];
        iterPara->startInner = 1;
		iterPara->endInner   = procData->xLength[0];
		iterPara->fixedValue = 1;
		break;
	case RIGHT:
        iterPara->startOuter = 1;
		iterPara->endOuter   = procData->xLength[2];
        iterPara->startInner = 1;
		iterPara->endInner   = procData->xLength[0];
		iterPara->fixedValue = procData->xLength[1];
		break;
	//---------------------------------------------
	//outer = Y, inner = X, Z fixed
	case TOP:
        iterPara->startOuter = 0;
		iterPara->endOuter   = procData->xLength[1]+1;
        iterPara->startInner = 1;
		iterPara->endInner   = procData->xLength[0];
		iterPara->fixedValue = procData->xLength[2];
		break;
	case BOTTOM:
        iterPara->startOuter = 0;
		iterPara->endOuter   = procData->xLength[1]+1;
        iterPara->startInner = 1;
		iterPara->endInner   = procData->xLength[0];
		iterPara->fixedValue = 1;
		break;

	//---------------------------------------------
	//outer = Z, inner = Y, X fixed
	case FRONT:
        iterPara->startOuter = 0;
		iterPara->endOuter   = procData->xLength[2]+1;
        iterPara->startInner = 0;
		iterPara->endInner   = procData->xLength[1]+1;
		iterPara->fixedValue = procData->xLength[0];
		break;
	case BACK:
        iterPara->startOuter = 0;
		iterPara->endOuter = procData->xLength[2]+1;
        iterPara->startInner = 0;
		iterPara->endInner = procData->xLength[1]+1;
		iterPara->fixedValue = 1;
		break;

	default:
		ERROR("Invalid direction occurred. This should never happen!");
	}
}

//Function to find indecies being extracted/injected
void p_assignIndices(int *direction, int *index) {
	switch (*direction) {
		case BOTTOM:
			// z = 0
			index[0] = 14; index[1] = 15; index[2] = 16; index[3] = 17; index[4] = 18;
			break;
		case TOP:
			// z = xlength[2]+1
			index[0] = 0; index[1] = 1; index[2] = 2; index[3] = 3; index[4] = 4;
			break;
		case BACK:
			// x = 0
			index[0] = 3; index[1] = 7; index[2] = 10; index[3] = 13; index[4] = 17;
			break;
		case FRONT:
			// x = xlength[0]+1
			index[0] = 1; index[1] = 5; index[2] = 8; index[3] =  11; index[4] = 15;
			break;
		case RIGHT:
			// y = xlength[1]+1
			index[0] = 0; index[1] = 5; index[2] = 6; index[3] = 7; index[4] = 14;
			break;
		case LEFT:
			// y = 0
			index[0] = 4; index[1] = 11; index[2] = 12; index[3] = 13; index[4] = 18;
			break;
	}
}

//Function to compute index in buffer given outer and inner coordinate
int p_computeBuffCellOffset(const int outer, const int inner, 
                            const int bufferLength[3][3], const int direction){


    //TODO: (TKS) Confirm that indecies are correct.
	//direction has valid integer values from 0 to 5
	switch (direction/2) { //integer division to get the type of face (see enum in LBDefinitions.h)
		case 0: // LEFT, RIGHT -> Y fixed
			//outer = Z, inner = X
			return bufferLength[direction/2][0]*(outer-1) + (inner-1);

		case 1: // TOP, BOTTOM -> Z fixed
			//outer = Y, inner = X
			return bufferLength[direction/2][0]*outer + (inner-1);

		case 2: // FRONT, BACK -> X fixed
			//outer = Z, inner = Y
			return bufferLength[direction/2][1]*outer + inner;

		default:
			ERROR("Invalid direction occured. This should not happen !!!");
			return -1;
	}
}
