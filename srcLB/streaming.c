#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
	/*
	 * For each FLUID field (this) copy distribution f_i from all neighboring
	 * cells of this.
	 */

	/* TODO: (DL) if there is time: optimize for cache efficiency.
	 * First: profile and see how this function adds to total workload (check if it's worth it)
	 */

	int zoffset, yzoffset, xyzoffset, current_cell_index; //current cell offset values
	int n_x, n_y, n_z, n_cell_index; 				  //neighbor offset values

	for(int z=1; z<=xlength; z++){
		zoffset = z*xlength;

		for(int y=1; y<=xlength; y++){
			yzoffset = xlength*(zoffset + y);

			for(int x=1; x<=xlength; x++){
				xyzoffset = yzoffset + x;

				/* TODO: (DL) if the indices are correct,
				* I think condition is always true... if so: delete if-statement
				*/
				// printf("Condition check %i \n", flagField[xyzoffset]);

				if( ! flagField[xyzoffset] ){ //true if FLUID cell
					current_cell_index = Q*xyzoffset;

					//loop through all neighbors and copy the respective distribution
					for(int i=0; i<Q; ++i){
						n_x = x+LATTICEVELOCITIES[i][0];
						n_y = y+LATTICEVELOCITIES[i][1];
						n_z = z+LATTICEVELOCITIES[i][2];
						n_cell_index = Q*(xlength*(n_z*xlength + n_y) + n_x);

						streamField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)];
					}
				}
			}
		}
	}
}
