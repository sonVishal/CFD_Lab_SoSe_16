#include "assert.h"
#include "computeCellValues.h"

/* Perform streaming and collision in one step */
void streamCollide(t_component *c, const t_procData * const procData) {

	// Define iteration indices
	int cellIdx;

	// Perform collision on all "inner" (FLUID) cells
	for (int z = 1; z <= procData->xLength[2] ; z++) {
		for (int y = 1; y <= procData->xLength[1]; y++) {
			for (int x = 1; x <= procData->xLength[0]; x++) {
				cellIdx = p_computeCellOffsetXYZ_Q(x, y, z, procData->xLength);
				for(int k = 0; k < numComp; ++k){
	                for (int i = 0; i < Q; ++i) {
						int nextCellIndex = p_computeCellOffsetXYZ_Q(x+LATTICEVELOCITIES[i][0],
							y+LATTICEVELOCITIES[i][1], z+LATTICEVELOCITIES[i][2], procData->xLength);

						c[k].field_new[cellIdx+Q-i-1] = c[k].field[nextCellIndex+Q-i-1]
                                                         - (c[k].field[nextCellIndex+Q-i-1]
                                                         - c[k].feq[nextCellIndex+Q-i-1])/c[k].tau;

						assert(c[k].field_new[cellIdx+Q-i-1] > 0.0);
	                }
				}

            }
        }
    }
}
