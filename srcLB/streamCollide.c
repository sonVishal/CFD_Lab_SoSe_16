#include "assert.h"
#include "computeCellValues.h"


void streamCollide(t_component *c, int xlength, double* feq, int* flagField){

	// Define iteration indices
	int cellIdx, fieldIdx;

	// Perform collision on all "inner" (FLUID) cells
	for (int z = 1; z <= xlength ; z++) {
		for (int y = 1; y <= xlength; y++) {
			for (int x = 1; x <= xlength; x++) {
				fieldIdx = p_computeCellOffsetXYZ(x, y, z, xlength);
				cellIdx = Q*fieldIdx;

				for(int k = 0; k < NUMCOMP; ++k){
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
}


void updateFeq(const int *xlength, const double*rho, double *velocity[3], double*feq){

	int cellIdx, fieldIdx;
    //TODO: fix components in main
    for(int k  = 1; k < NUMCOMP; ++k){
        for (int z = 1; z <= *xlength ; z++) {
            for (int y = 1; y <= *xlength; y++) {
                for (int x = 1; x <= *xlength; x++) {

                    fieldIdx = p_computeCellOffsetXYZ(x, y, z, *xlength);
                    cellIdx = Q*fieldIdx;

                    computeFeq(&rho[fieldIdx], velocity[fieldIdx], &feq[cellIdx]);


                }
            }
        }
    }
}
