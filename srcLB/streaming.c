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

    //TODO: (TKS) Found the whole offset notation a bit confusing.
    //            make a variable called index in the inner loop. That feels
    //            more descriptive. index gives the index of the first 
    //            direction in a cell at position (x,y,z)
    //
    //            Suggested name changes:
    //            cell_xyzoffset --> index or current_cell_index
    //            n_xyzoffset --> n_index or n_cell_index
    //
    //

	int zoffset, yzoffset, xyzoffset, cell_xyzoffset; //current cell offset values
	int n_x, n_y, n_z, n_xyzoffset; 				  //neighbor offset values

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
					cell_xyzoffset = Q*xyzoffset;

					//loop through all neighbors and copy the respective distribution
					for(int i=0; i<Q; ++i){
						n_x = x+LATTICEVELOCITIES[i][0];
						n_y = y+LATTICEVELOCITIES[i][1];
						n_z = z+LATTICEVELOCITIES[i][2];
						n_xyzoffset = Q*(xlength*(n_z*xlength + n_y) + n_x);

						streamField[cell_xyzoffset + i] = collideField[n_xyzoffset + (Q-i-1)];
					}
				}
			}
		}
	}
}
