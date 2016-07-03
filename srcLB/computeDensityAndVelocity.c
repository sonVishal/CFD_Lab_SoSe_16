#include "LBDefinitions.h"
#include "computeCellValues.h"

//void computeVelocity(double *currentCell, double *velocity, double *density){

    //velocity[0] = 0;
    //velocity[1] = 0;
    //velocity[2] = 0;

	//for (int i = 0; i < Q; i++) {
		   //velocity[0] += currentCell[i]*LATTICEVELOCITIES[i][0];
		   //velocity[1] += currentCell[i]*LATTICEVELOCITIES[i][1];
		   //velocity[2] += currentCell[i]*LATTICEVELOCITIES[i][2];
	//}

    //velocity[0] /= (*density);
    //velocity[1] /= (*density);
    //velocity[2] /= (*density);
//}

//void computeNumDensity(double *currentCell, double *density){
    //// Density is the sum of the distributions in the current lattice
    //int i;
    //*density = 0.0;
	//for (i = 0; i < Q; i++) {
		 //*density += currentCell[i];
	//}
//}

void computeDensityAndVelocity(t_component *c, int xlength){

    for (int z = 1; z <= xlength ; z++) {
        for (int y = 1; y <= xlength; y++) {
            for (int x = 1; x <= xlength; x++) {

                for(int k = 0; k < NUMCOMP; ++k){
                    int fieldIdx = p_computeCellOffsetXYZ(x, y, z, xlength);
                    int cellIdx = Q*fieldIdx;
                    double *currentCell = &c[k].collideField[cellIdx];

                    //TODO: this has to be done component wise..., actually velocity, density and force have to be saved in the component struct
                    for(int k = 0; k < NUMCOMP; ++k){
                        computeNumDensity(currentCell, &c[k].rho[fieldIdx]);
                        computeVelocity(currentCell, c[k].velocity[fieldIdx], &c[k].rho[fieldIdx]);
                    }

                }
            }
        }
    }
}
