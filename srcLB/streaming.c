#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
	/*
	 * For each FLUID field (this) copy distribution f_i from all neighboring
	 * cells of this.
	 */

	int xlen_2 = xlength+2;
	int xlen_2sq = xlen_2*xlen_2;
	int totalSize = xlen_2*xlen_2sq;

	int i;

	int nextCellIndex, currentCellIndex;

	for (i = 0; i < totalSize; i++) {
		if (flagField[i] == 0) {
			currentCellIndex = Q*i;

			// Loop unroll
			nextCellIndex = Q*(i+xlen_2sq+xlen_2);
			streamField[currentCellIndex] = collideField[nextCellIndex];

			nextCellIndex = Q*(i+xlen_2sq+1);
			streamField[currentCellIndex+1] = collideField[nextCellIndex+1];

			nextCellIndex = Q*(i+xlen_2sq);
			streamField[currentCellIndex+2] = collideField[nextCellIndex+2];

			nextCellIndex = Q*(i+xlen_2sq-1);
			streamField[currentCellIndex+3] = collideField[nextCellIndex+3];

			nextCellIndex = Q*(i+xlen_2sq-xlen_2);
			streamField[currentCellIndex+4] = collideField[nextCellIndex+4];

			nextCellIndex = Q*(i+xlen_2+1);
			streamField[currentCellIndex+5] = collideField[nextCellIndex+5];

			nextCellIndex = Q*(i+xlen_2);
			streamField[currentCellIndex+6] = collideField[nextCellIndex+6];

			nextCellIndex = Q*(i+xlen_2-1);
			streamField[currentCellIndex+7] = collideField[nextCellIndex+7];

			nextCellIndex = Q*(i+1);
			streamField[currentCellIndex+8] = collideField[nextCellIndex+8];

			streamField[currentCellIndex+9] = collideField[currentCellIndex+9];

			nextCellIndex = Q*(i-1);
			streamField[currentCellIndex+10] = collideField[nextCellIndex+10];

			nextCellIndex = Q*(i+1-xlen_2);
			streamField[currentCellIndex+11] = collideField[nextCellIndex+11];

			nextCellIndex = Q*(i-xlen_2);
			streamField[currentCellIndex+12] = collideField[nextCellIndex+12];

			nextCellIndex = Q*(i-1-xlen_2);
			streamField[currentCellIndex+13] = collideField[nextCellIndex+13];

			nextCellIndex = Q*(i+xlen_2-xlen_2sq);
			streamField[currentCellIndex+14] = collideField[nextCellIndex+14];

			nextCellIndex = Q*(i+1-xlen_2sq);
			streamField[currentCellIndex+15] = collideField[nextCellIndex+15];

			nextCellIndex = Q*(i-xlen_2sq);
			streamField[currentCellIndex+16] = collideField[nextCellIndex+16];

			nextCellIndex = Q*(i-1-xlen_2sq);
			streamField[currentCellIndex+17] = collideField[nextCellIndex+17];

			nextCellIndex = Q*(i-xlen_2-xlen_2sq);
			streamField[currentCellIndex+18] = collideField[nextCellIndex+18];
		}
	}
}
