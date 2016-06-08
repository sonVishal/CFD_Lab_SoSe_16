
#include "helperMPI.h"

//TODO: (TKS) 
//Reduce bufferSize[6]--> bufferSize[3]?
//          * Can do two iterations of the loop at the same time.
//          * Or complete loop unroll
//
//At what level to iterate plane?
//          * inside each function
//          * make a double for in communicate and use inner/outer formulation.
//                  * Extract/Inject need only get indecies of cillideField + buffer.
//

void p_assignIndices(int *direction, int *index);

void communicate(double** sendBuffer, double**readBuffer, double* collideField, int* xlength, int* bufferSize){
    //Run extract, swap, inject for all sides and cells.
    //
    int index[5];
    for (int direction = LEFT; direction < BACK; ++direction) {
        //TODO: (TKS) iterate over all planes and do the following:
        p_assignIndices(&direction, index);
        extract(sendBuffer, collideField, xlength, bufferSize, direction);
        swap(sendBuffer, readBuffer, bufferSize, direction);
        inject(readBuffer, collideField, xlength, bufferSize, direction);
    }

}

//Copy distributions needed to sendbuffer.
void extract( double** sendBuffer, double* collideField, int* xlength, int* bufferSize, int direction){
}

//Send distributions and wait to receive.
void swap(double** sendBuffer, double**readBuffer, int* bufferSize, int direction){

}

//Copy read buffer to ghost layer
void inject(double** readBuffer, double* collideField, int* xlength, int* bufferSize, int direction){
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
