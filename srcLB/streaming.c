#include "streaming.h"


/* Another (obsolete) version where we iterated over FLUID cells with 3 loops (x,y,z) turned out to
 * to be slower than the current version.
 *
 * Here we just iterate over flagField and depending on the value (if FLUID cell) carry out the streaming.
 * Furthermore, the loop is unrolled, which gained better performance.
 */
void doStreaming(t_component *c,int *flagField,int xlength){
	/*
	 * Each FLUID cell (this) copies the distribution f_i-values of all neighboring
	 * cells of this.
	 */

	int nextCellIndex, currentCellIndex;

	int x,y,z;
	for (int k = 0; k < NUMCOMP; k++) {
		for (z = 1; z <= xlength; z++) {
			for (y = 1; y <= xlength; y++) {
				for (x = 1; x <= xlength; x++) {
					currentCellIndex = p_computeCellOffsetXYZ_Q(x, y, z, xlength);
					for (int i = 0; i < Q; i++) {
						int nextPoint[3];
						nextPoint[0] = (x+LATTICEVELOCITIES[i][0] <= xlength) ? x+LATTICEVELOCITIES[i][0]:1;
			            nextPoint[0] = (x+LATTICEVELOCITIES[i][0] >= 0) ? x+LATTICEVELOCITIES[i][0]:xlength;
			            nextPoint[1] = (y+LATTICEVELOCITIES[i][1] <= xlength) ? y+LATTICEVELOCITIES[i][1]:1;
			            nextPoint[1] = (y+LATTICEVELOCITIES[i][1] >= 0) ? y+LATTICEVELOCITIES[i][1]:xlength;
			            nextPoint[2] = (z+LATTICEVELOCITIES[i][2] <= xlength) ? z+LATTICEVELOCITIES[i][2]:1;
			            nextPoint[2] = (z+LATTICEVELOCITIES[i][2] >= 0) ? z+LATTICEVELOCITIES[i][2]:xlength;
						nextCellIndex = p_computeCellOffsetXYZ_Q(nextPoint[0], nextPoint[1], nextPoint[2], xlength);
						c[k].streamField[currentCellIndex+i] = c[k].collideField[nextCellIndex+i];
					}
				}
			}
		}
	}

	// for (int k = 0; k < NUMCOMP; k++) {
	// // Inner loop unrolling
	// // Faster even with -O3
	// // Up to 2 seconds with xlength = 30
	// 	for (i = 0; i < totalSize; i++) {
	// 		if (flagField[i] == 0) {
	// 			currentCellIndex = Q*i;
	//
	// 			// Loop unroll
	// 			//current cell to neighbor cell in direction LATTICEVELOCITIES[18] = {0,1,1}
	// 			//distribution j = 0
	// 			nextCellIndex = Q*(i+xlen_2sq+xlen_2);
	// 			c[k].streamField[currentCellIndex] = c[k].collideField[nextCellIndex];
	//
	// 			//current cell to neighbor cell in direction LATTICEVELOCITIES[17] = {1,0,1}
	// 			//distribution j = 1
	// 			nextCellIndex = Q*(i+xlen_2sq+1);
	// 			c[k].streamField[currentCellIndex+1] = c[k].collideField[nextCellIndex+1];
	//
	// 			//current cell to neighbor cell in direction LATTICEVELOCITIES[16] = {0,0,1}
	// 			//distribution j = 2
	// 			nextCellIndex = Q*(i+xlen_2sq);
	// 			c[k].streamField[currentCellIndex+2] = c[k].collideField[nextCellIndex+2];
	//
	// 			//current cell to neighbor cell in direction LATTICEVELOCITIES[Q-j-1]
	// 			//distribution j = 3
	// 			nextCellIndex = Q*(i+xlen_2sq-1);
	// 			c[k].streamField[currentCellIndex+3] = c[k].collideField[nextCellIndex+3];
	//
	// 			//etc...
	// 			nextCellIndex = Q*(i+xlen_2sq-xlen_2);
	// 			c[k].streamField[currentCellIndex+4] = c[k].collideField[nextCellIndex+4];
	//
	// 			nextCellIndex = Q*(i+xlen_2+1);
	// 			c[k].streamField[currentCellIndex+5] = c[k].collideField[nextCellIndex+5];
	//
	// 			nextCellIndex = Q*(i+xlen_2);
	// 			c[k].streamField[currentCellIndex+6] = c[k].collideField[nextCellIndex+6];
	//
	// 			nextCellIndex = Q*(i+xlen_2-1);
	// 			c[k].streamField[currentCellIndex+7] = c[k].collideField[nextCellIndex+7];
	//
	// 			nextCellIndex = Q*(i+1);
	// 			c[k].streamField[currentCellIndex+8] = c[k].collideField[nextCellIndex+8];
	//
	// 			c[k].streamField[currentCellIndex+9] = c[k].collideField[currentCellIndex+9];
	//
	// 			nextCellIndex = Q*(i-1);
	// 			c[k].streamField[currentCellIndex+10] = c[k].collideField[nextCellIndex+10];
	//
	// 			nextCellIndex = Q*(i+1-xlen_2);
	// 			c[k].streamField[currentCellIndex+11] = c[k].collideField[nextCellIndex+11];
	//
	// 			nextCellIndex = Q*(i-xlen_2);
	// 			c[k].streamField[currentCellIndex+12] = c[k].collideField[nextCellIndex+12];
	//
	// 			nextCellIndex = Q*(i-1-xlen_2);
	// 			c[k].streamField[currentCellIndex+13] = c[k].collideField[nextCellIndex+13];
	//
	// 			nextCellIndex = Q*(i+xlen_2-xlen_2sq);
	// 			c[k].streamField[currentCellIndex+14] = c[k].collideField[nextCellIndex+14];
	//
	// 			nextCellIndex = Q*(i+1-xlen_2sq);
	// 			c[k].streamField[currentCellIndex+15] = c[k].collideField[nextCellIndex+15];
	//
	// 			nextCellIndex = Q*(i-xlen_2sq);
	// 			c[k].streamField[currentCellIndex+16] = c[k].collideField[nextCellIndex+16];
	//
	// 			nextCellIndex = Q*(i-1-xlen_2sq);
	// 			c[k].streamField[currentCellIndex+17] = c[k].collideField[nextCellIndex+17];
	//
	// 			nextCellIndex = Q*(i-xlen_2-xlen_2sq);
	// 			c[k].streamField[currentCellIndex+18] = c[k].collideField[nextCellIndex+18];
	// 		}
	// 	}
	// }
}
