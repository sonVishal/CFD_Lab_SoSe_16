#include "LBDefinitions.h"



void streamCollide(int xlength, double* streamField, double* collideField, double* feq, int* flagField){

	// Define iteration indices
	int cellIdx, fieldIdx;

	// Perform collision on all "inner" (FLUID) cells
	for (int z = 1; z <= xlength ; z++) {
		for (int y = 1; y <= xlength; y++) {
			for (int x = 1; x <= xlength; x++) {

				// Get the index of the first distribution
				// in the current cell
				fieldIdx = p_computeCellOffsetXYZ(x, y, z, xlength);
				cellIdx = Q*fieldIdx;

                for (int i = 0; i < Q; ++i) {
					//TODO: Check if the order is correct

					// f_new[f_index_end] = f[f_index_beg]
					//                    - (f[f_index_beg] - f_eq[f_index_beg])
					//                    / tau;
					double tau = 1.0;
			        collideField[cellIdx+Q-i-1] = streamField[cellIdx+i] - (streamField[cellIdx+i] - feq[cellIdx+i]) / tau;
                }
            }
        }
    }
}
