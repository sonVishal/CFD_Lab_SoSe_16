#include "assert.h"
#include "computeCellValues.h"


void streamCollide(t_component *c, int xlength, int* flagField){

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
						int nextCellIndex = x+LATTICEVELOCITIES[i][0] + (y+LATTICEVELOCITIES[i][1])*(xlength+2) +
											(z+LATTICEVELOCITIES[i][2])*(xlength+2)*(xlength+2);
						// assert(c[k].tau == 1.0);
						c[k].collideField[cellIdx+Q-i-1] = c[k].streamField[nextCellIndex*Q+Q-i-1] - (c[k].streamField[nextCellIndex*Q+Q-i-1] - c[k].feq[nextCellIndex*Q+Q-i-1])/c[k].tau;
						assert(c[k].collideField[cellIdx+Q-i-1] > 0.0);
	                }
				}

            }
        }
    }
}


void updateFeq(t_component *c, const int *xlength){

	int cellIdx, fieldIdx;
    //TODO: fix components in main
    for(int k  = 0; k < NUMCOMP; ++k){
        for (int z = 1; z <= *xlength ; z++) {
            for (int y = 1; y <= *xlength; y++) {
                for (int x = 1; x <= *xlength; x++) {

                    fieldIdx = p_computeCellOffsetXYZ(x, y, z, *xlength);
                    cellIdx = Q*fieldIdx;
                    double ux = c[k].velocity[fieldIdx][0];
                    double uy = c[k].velocity[fieldIdx][1];
                    double uz = c[k].velocity[fieldIdx][2];

                    double udotu = ux*ux + uy*uy + uz*uz;

                    //computeFeq(&c[k].rho[fieldIdx], c[k].velocity[fieldIdx], &c[k].feq[cellIdx]);
					for (int i = 0; i < Q; i++) {
                        double edotu =  LATTICEVELOCITIES[i][0]*ux +
                                        LATTICEVELOCITIES[i][1]*uy +
                                        LATTICEVELOCITIES[i][2]*uz;

                        c[k].feq[cellIdx + i] = LATTICEWEIGHTS[i] * c[k].rho[fieldIdx]
                                                      * (1 + 3.*edotu + 4.5*edotu*edotu - 1.5*udotu);
						assert(c[k].feq[cellIdx+i] > 0.0);
					}

                }
            }
        }
    }
}
