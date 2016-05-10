#include "streaming.h"
#include "LBDefinitions.h"
#include <stdio.h>

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
	/*
	 * For each FLUID field (this) copy distribution f_i from all neighboring
	 * cells of this.
	 */

	int zoffset, yzoffset, xyzoffset;  // Temporary variables to save computation
	int current_cell_index;            // Position of first direction entry of the cell
	// positioned at (x,y,z)

	int n_x, n_y, n_z, n_cell_index;  // Neighbor offset values

	for(int z=1; z<=xlength; z++){
		zoffset = z*(xlength+2);

		for(int y=1; y<=xlength; y++){
			yzoffset = (xlength+2)*(zoffset + y);

			for(int x=1; x<=xlength; x++){
				xyzoffset = yzoffset + x;

				//Code snippet in case of debugging. We dont check for FLUID condition because it is
				//always true (iterating only FLUID cells) by setting the indices correctly.
//				if(flagField[xyzoffset]){
//					printf("Condition check %i \n",flagField[xyzoffset] );
//				}

				current_cell_index = Q*xyzoffset;

				// Loop through all neighbors and copy their respective
				// distributions to the streamField.
				for(int i=0; i<Q; ++i){
					n_x = x-LATTICEVELOCITIES[i][0];
					n_y = y-LATTICEVELOCITIES[i][1];
					n_z = z-LATTICEVELOCITIES[i][2];
					n_cell_index = Q*((xlength+2)*(n_z*(xlength+2) + n_y) + n_x);

					streamField[current_cell_index + i] = collideField[n_cell_index + i];
				}

			}
		}
	}
}
