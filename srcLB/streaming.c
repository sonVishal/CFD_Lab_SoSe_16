#include <stdio.h>

#include "streaming.h"
#include "LBDefinitions.h"



/* Another (obsolete) version where we iterated over FLUID cells with 3 loops (x,y,z) turned out to
 * to be slower than the current version.
 *
 * Here we just iterate over flagField and depending on the value (if FLUID cell) carry out the streaming.
 * Furthermore, the loop is unrolled, which gained better performance.
 */
void doStreaming(double *collideField, double *streamField,int *flagField,int *xlength){
	/*
	 * Each FLUID cell (this) copies the distribution f_i-values of all neighboring
	 * cells of this.
	 */

	/*TODO:(DL) Check if this is correct after finishing with initialization
	* adapted to xlength[3] */

    int xlen2 = xlength[0]+2;
	int ylen2 = xlength[1]+2;
    int xylen2 = xlen2*ylen2;

	int zlen2 = xlength[1]+2;
	int totalSize = xlen2*ylen2*zlen2;
	int i;

	int nextCellIndex, currentCellIndex;

	// Semantics for the loop unrolling
	// int i, j;
	// for (i = 0; i < totalSize; i++) {
	// 	if (flagField[i] == 0) {
	// 		currentCellIndex = Q*i;
	// 		for (j = 0; j < Q; j++) {
	// 			nextCellIndex = Q*(i-LATTICEVELOCITIES[j][0]-
	// 				xlen_2*LATTICEVELOCITIES[j][1]-
	// 				xlen_2sq*LATTICEVELOCITIES[j][2]);
	// 			streamField[currentCellIndex+j] = collideField[nextCellIndex+j];
	// 		}
	// 	}
	// }

	// Inner loop unrolling
	// Faster even with -O3
	// Up to 2 seconds with xlength = 30

	for (i = 0; i < totalSize; i++) {

		if (flagField[i] == 0) {
			currentCellIndex = Q*i;

			// Loop unroll
			//current cell to neighbor cell in direction LATTICEVELOCITIES[18] = {0,1,1}
			//distribution j = 0

			nextCellIndex = Q*(i+xylen2+xlen2);
			streamField[currentCellIndex] = collideField[nextCellIndex];

			//current cell to neighbor cell in direction LATTICEVELOCITIES[17] = {1,0,1}
			//distribution j = 1
			nextCellIndex = Q*(i+xylen2+1);
			streamField[currentCellIndex+1] = collideField[nextCellIndex+1];

			//current cell to neighbor cell in direction LATTICEVELOCITIES[16] = {0,0,1}
			//distribution j = 2
			nextCellIndex = Q*(i+xylen2);
			streamField[currentCellIndex+2] = collideField[nextCellIndex+2];

			//current cell to neighbor cell in direction LATTICEVELOCITIES[Q-j-1]
			//distribution j = 3
			nextCellIndex = Q*(i+xylen2-1);
			streamField[currentCellIndex+3] = collideField[nextCellIndex+3];

			//etc...
			nextCellIndex = Q*(i+xylen2-xlen2);
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

			nextCellIndex = Q*(i+xlen2-xylen2);
			streamField[currentCellIndex+14] = collideField[nextCellIndex+14];

			nextCellIndex = Q*(i+1-xylen2);
			streamField[currentCellIndex+15] = collideField[nextCellIndex+15];

			nextCellIndex = Q*(i-xylen2);
			streamField[currentCellIndex+16] = collideField[nextCellIndex+16];

			nextCellIndex = Q*(i-1-xylen2);
			streamField[currentCellIndex+17] = collideField[nextCellIndex+17];

			nextCellIndex = Q*(i-xlen2-xylen2);
			streamField[currentCellIndex+18] = collideField[nextCellIndex+18];
		}
	}
}
