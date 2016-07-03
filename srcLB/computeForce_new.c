
#include "LBDefinitions.h"
#include "assert.h"


void computeCellForce(const int currentCellIndex, const int currentCompIndex,
    const t_component *const c, const int * const flagField,
    double const*const G, int xlength, double forces[3]) {

    // Important: The currentCellIndex is without multiplication with Q
    int xlen2 = xlength+2;
    int xlen2sq = xlen2*xlen2;

    double numDensity;
    forces[0] = 0.0; forces[1] = 0.0; forces[2] = 0.0;

    for (int m = 0; m < NUMCOMP; m++) {
        for (int i = 0; i < Q; i++) {

            // Go to the next cell index in the direction of lattice velocities
            int nextCellIndex = currentCellIndex+LATTICEVELOCITIES[i][0]
                                + xlen2*LATTICEVELOCITIES[i][1]
                                + xlen2sq*LATTICEVELOCITIES[i][2]; //index of cell in direction i

            // int nextIndex = Q*nextCellIndex; //index of number density in direction i
            // numDensity = c[m].collideField[nextIndex]; //number density in direction i
            // Compute the number density for component "m" at lattice site "nextIndex"

            // computeNumDensity(&c[m].collideField[nextIndex], &numDensity);
            numDensity = c[m].rho[nextCellIndex];
            // if(flagField[nextCellIndex] == FLUID){
            // }else{
            //     numDensity = c[m].collideField[nextIndex + 9];
            // }

            double G_cur;
            if(LATTICEWEIGHTS[i] == w2){ // 1/18
                G_cur = G[m];
                assert(G_cur == -0.27);
            }else if(LATTICEWEIGHTS[i] == w3){ // 1/36
                G_cur = G[m]/2;
                assert(G_cur == -0.27/2);
            }else{
                assert(i == 9);
                G_cur = 0;
            }

            //Shan&Doolen eq. 4 (PDF page 5)
            forces[0] += G_cur * psiFctPointer[c[m].psiFctCode](numDensity) * LATTICEVELOCITIES[i][0];
            forces[1] += G_cur * psiFctPointer[c[m].psiFctCode](numDensity) * LATTICEVELOCITIES[i][1];
            forces[2] += G_cur * psiFctPointer[c[m].psiFctCode](numDensity) * LATTICEVELOCITIES[i][2];
            assert(psiFctPointer[c[m].psiFctCode](numDensity) == psi0(numDensity));
         }
    }

    // numDensity = c[n].collideField[currentCellIndex];
    // computeNumDensity(&c[currentCompIndex].collideField[currentCellIndex], &numDensity);
    numDensity = c[currentCompIndex].rho[currentCellIndex];

    forces[0] *= -psiFctPointer[c[currentCompIndex].psiFctCode](numDensity);
    forces[1] *= -psiFctPointer[c[currentCompIndex].psiFctCode](numDensity);
    forces[2] *= -psiFctPointer[c[currentCompIndex].psiFctCode](numDensity);
    assert(currentCompIndex == 0);
    assert(psiFctPointer[c[currentCompIndex].psiFctCode](numDensity) == psi0(numDensity));
}

void computeForce_new(t_component *c, int xlength, int *flagField, double **G){

    for (int z = 1; z <= xlength ; z++) {
        for (int y = 1; y <= xlength; y++) {
            for (int x = 1; x <= xlength; x++) {
                int fieldIdx = p_computeCellOffsetXYZ(x, y, z, xlength);
                // int cellIdx = Q*fieldIdx;
                for(int k = 0; k < NUMCOMP; ++k){
                    computeCellForce(fieldIdx, k, c, flagField, G[k], xlength, c[k].force[fieldIdx]);
                }

            }
        }
    }

}
