#include "LBDefinitions.h"



void streamCollide(t_component *c, int xlength, double* feq, int* flagField){

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
					// f_new[f_index_end] = f[f_index_beg]
					//                    - (f[f_index_beg] - f_eq[f_index_beg])
					//                    / tau;
					//TODO: Check if this is correct:
					c->collideField[cellIdx+Q-i-1] = c->streamField[cellIdx+i] - (c->streamField[cellIdx+i] - feq[cellIdx+i])/c->tau;
                }
            }
        }
    }
}
