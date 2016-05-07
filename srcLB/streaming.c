#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){

	/*
	 * For each FLUID field (this) copy distribution f_i from all neighboring
	 * cells.
	 */

	for(int z=1; z<=xlength; z++){
		int zoffset = z*xlength;

		for(int y=1; y<=xlength; y++){
			int yzoffset = xlength*(zoffset + y);

			for(int x=1; x<=xlength; x++){
				int xyzoffset = yzoffset + x;

				// TODO: (DL) if the indices are correct, I think condition is always true...
				// printf("Condition check \i", flagField[xyzoffset]);

				if( ! flagField[xyzoffset] ){ //true if FLUID cell
					int dist_xyzoffset = 19*xyzoffset;

					//loop through all neighbors and copy the respective distribution
					for(int i=0; i<19; ++i){
						streamField[dist_xyzoffset + i] = collideField[dist_xyzoffset + (19-i)];
					}
				}
			}
		}
	}
}
