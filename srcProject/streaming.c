#include "streaming.h"

//TODO:(TKS) Make wrapper for several components
void streamComponents(t_component* c, int *flagField, int *xlength){
    for (int i = 0; i < numComp; ++i) {
        doStreaming(c[i].streamField, c[i].collideField, flagField, xlength);
    }
}

/* Another (obsolete) version where we iterated over FLUID cells with 3 loops (x,y,z) turned out to
 * to be slower than the current version.
 *
 * Here we just iterate over flagField and depending on the value (if FLUID cell) carry out the streaming.
 * Furthermore, the loop is unrolled, which gained better performance.
 */
void doStreaming(double *streamField, double *collideField, int *flagField,int *xlength){
	/*
	 * Each FLUID cell (this) copies the distribution f_i-values of all neighboring
	 * cells of this.
	 */

    int xlen2 = xlength[0]+2;
	int ylen2 = xlength[1]+2;
    int xylen = xlen2*ylen2;

	int zlen2 = xlength[2]+2;
	int totalSize = xlen2*ylen2*zlen2;
	int i;

	int nextCellIndex, currentCellIndex;

	// Semantics for the loop unrolling
	// int j, nextFlagIndex;
	// for (i = 0; i < totalSize; i++) {
	// 	if (flagField[i] == FLUID) {
	// 		currentCellIndex = Q*i;
	// 		for (j = 0; j < Q; j++) {
	// 			nextFlagIndex = i-LATTICEVELOCITIES[j][0]-
	// 							xlen2*LATTICEVELOCITIES[j][1]-
	// 							xlen2*ylen2*LATTICEVELOCITIES[j][2];
	// 				nextCellIndex = Q*nextFlagIndex;
	// 				streamField[currentCellIndex+j] = collideField[nextCellIndex+j];
	// 		}
	// 	}
	// }

	// Inner loop unrolling
	// Faster even with -O3
	// Up to 2 seconds with xlength = 30

	for (i = 0; i < totalSize; i++) {
		if (flagField[i] == FLUID) {
			currentCellIndex = Q*i;

			// Loop unroll
			//current cell to neighbor cell in direction LATTICEVELOCITIES[18] = {0,1,1}
			//distribution j = 0

			nextCellIndex = Q*(i+xylen+xlen2);
            streamField[currentCellIndex] = collideField[nextCellIndex];

			//current cell to neighbor cell in direction LATTICEVELOCITIES[17] = {1,0,1}
			//distribution j = 1
			nextCellIndex = Q*(i+xylen+1);
            streamField[currentCellIndex+1] = collideField[nextCellIndex+1];

			//current cell to neighbor cell in direction LATTICEVELOCITIES[16] = {0,0,1}
			//distribution j = 2
			nextCellIndex = Q*(i+xylen);
            streamField[currentCellIndex+2] = collideField[nextCellIndex+2];

			//current cell to neighbor cell in direction LATTICEVELOCITIES[Q-j-1]
			//distribution j = 3
			nextCellIndex = Q*(i+xylen-1);
            streamField[currentCellIndex+3] = collideField[nextCellIndex+3];

			//etc...
			nextCellIndex = Q*(i+xylen-xlen2);
            streamField[currentCellIndex+4] = collideField[nextCellIndex+4];

			nextCellIndex = Q*(i+xlen2+1);
            streamField[currentCellIndex+5] = collideField[nextCellIndex+5];

			nextCellIndex = Q*(i+xlen2);
            streamField[currentCellIndex+6] = collideField[nextCellIndex+6];

			nextCellIndex = Q*(i+xlen2-1);
            streamField[currentCellIndex+7] = collideField[nextCellIndex+7];

			nextCellIndex = Q*(i+1);
            streamField[currentCellIndex+8] = collideField[nextCellIndex+8];

			streamField[currentCellIndex+9] = collideField[currentCellIndex+9];

			nextCellIndex = Q*(i-1);
            streamField[currentCellIndex+10] = collideField[nextCellIndex+10];

			nextCellIndex = Q*(i+1-xlen2);
            streamField[currentCellIndex+11] = collideField[nextCellIndex+11];

			nextCellIndex = Q*(i-xlen2);
            streamField[currentCellIndex+12] = collideField[nextCellIndex+12];

			nextCellIndex = Q*(i-1-xlen2);
            streamField[currentCellIndex+13] = collideField[nextCellIndex+13];

			nextCellIndex = Q*(i+xlen2-xylen);
            streamField[currentCellIndex+14] = collideField[nextCellIndex+14];

			nextCellIndex = Q*(i+1-xylen);
            streamField[currentCellIndex+15] = collideField[nextCellIndex+15];

			nextCellIndex = Q*(i-xylen);
            streamField[currentCellIndex+16] = collideField[nextCellIndex+16];

			nextCellIndex = Q*(i-1-xylen);
            streamField[currentCellIndex+17] = collideField[nextCellIndex+17];

			nextCellIndex = Q*(i-xlen2-xylen);
            streamField[currentCellIndex+18] = collideField[nextCellIndex+18];
		}
	}
}
