
#include "helperMPI.h"
#include "boundary.h"
#include "helper.h"

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
    int startInner, endInner;
    int startOuter, endOuter;
    int fixedValue;

    for (int direction = LEFT; direction < BACK; ++direction) {
        //TODO: (TKS) iterate over all planes and do the following:
        //          * Using inner/outer formulation.
        //          * Functions operate on single cells.
        //              - Need functions to take the whole extraction etc?

        p_setCommIterationParameters(&startOuter, &endOuter,&startInner, &endInner, &fixedValue, *procData, direction);
        p_assignIndices(&direction, index);

        //TODO:(TKS) Should move for loops into the functions. Want to do whole extraction before doing swap etc.
	    //k - corresponds to the 'outer' value when computing the offset
        for(int k = startOuter; k <= endOuter; ++k){
            //j - corresponds to the 'inner' value
            for(int j = startInner; j <= endInner; ++j){

                int currentCellIndex = Q*p_computeCellOffset(k, j, fixedValue, (*procData).xLength, direction);

                //TODO: (TKS) Decide whether to take in procData of introduce temp variables, xlength, bufferSize.
                extract(sendBuffer, collideField, procData, direction);
                swap(sendBuffer, readBuffer, procData, direction);
                inject(readBuffer, collideField, procData, direction);
            }
        }
    }

}

//Copy distributions needed to sendbuffer.
void extract( double** sendBuffer, double* collideField, const t_procData *procData, int direction){
}

//Send distributions and wait to receive.
void swap(double** sendBuffer, double**readBuffer, const t_procData *procData, int direction){

}

//Copy read buffer to ghost layer
void inject(double** readBuffer, double* collideField, const t_procData *procData, int direction){
}

//Function to assign iteration parameters for communication.

void p_setCommIterationParameters(int *startOuter, int *endOuter, int *startInner, int *endInner, 
                                  int *fixedValue, const t_procData procData, const int direction){

    //TODO: (TKS) Confirm that indecies are correct.
	switch(direction){

	//---------------------------------------------
	//outer = Z, inner = X, Y fixed
    //only iterate over inner domain of plane (FLUID cells)
	case LEFT:
        *startOuter = 1;
		*endOuter   = procData.xLength[2];
        *startInner = 1;
		*endInner   = procData.xLength[0];
		*fixedValue = 1;
		break;
	case RIGHT:
        *startOuter = 1;
		*endOuter   = procData.xLength[2];
        *startInner = 1;
		*endInner   = procData.xLength[0];
		*fixedValue = procData.xLength[1];
		break;
	//---------------------------------------------
	//outer = Y, inner = X, Z fixed
	case TOP:
        *startOuter = 0;
		*endOuter   = procData.xLength[1]+1;
        *startInner = 1;
		*endInner   = procData.xLength[0];
		*fixedValue = procData.xLength[2];
		break;
	case BOTTOM:
        *startOuter = 0;
		*endOuter   = procData.xLength[1]+1;
        *startInner = 1;
		*endInner   = procData.xLength[0];
		*fixedValue = 1;
		break;

	//---------------------------------------------
	//outer = Z, inner = Y, X fixed
	case FRONT:
        *startOuter = 0;
		*endOuter   = procData.xLength[2]+1;
        *startInner = 0;
		*endInner   = procData.xLength[1]+1;
		*fixedValue = procData.xLength[0];
		break;
	case BACK:
        *startOuter = 0;
		*endOuter = procData.xLength[2]+1;
        *startInner = 0;
		*endInner = procData.xLength[1]+1;
		*fixedValue = 1;
		break;

	default:
		ERROR("Invalid wallIdx occurred. This should never happen!");
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
