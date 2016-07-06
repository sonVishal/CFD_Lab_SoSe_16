#include "computeCellValues.h"

/** computes the density from the particle distribution functions stored at
 *  currentCell. currentCell thus denotes the address of the first particle
 *  distribution function of the respective cell.
 *  The result is stored in density.
 */
void computeNumDensity(const double *const currentCell, double *density){
    // Density is the sum of the distributions in the current lattice
    int i;
    *density = 0.0;
    for (i = 0; i < Q; i++) {
        *density += currentCell[i];
    }
}

/** computes the velocity within currentCell and stores the result in velocity */
void computeVelocity(const double * const currentCell, const double * const density, double *velocity){

    // Velocity is the momentum divided by the density
    // Momentum is the sum of the product of lattice velocity with distribution

    // Semantics of the unrolled loop
    // velocity[0] = 0.0;
    // velocity[1] = 0.0;
    // velocity[2] = 0.0;
    // int i;
    // for (i = 0; i < Q; i++) {
    // 	   velocity[0] += currentCell[i]*LATTICEVELOCITIES[i][0];
    // 	   velocity[1] += currentCell[i]*LATTICEVELOCITIES[i][1];
    // 	   velocity[2] += currentCell[i]*LATTICEVELOCITIES[i][2];
    // }

    // Unroll the loop for cache efficient access
    // Improved speed even with -O3 flag
    if ((*density) != 0.0) {
        // Blocks of 4 for efficient cache access
        // 0 to 3
        velocity[0] = -currentCell[1]+currentCell[3];
        velocity[1] = -currentCell[0];
        velocity[2] = -(currentCell[0]+currentCell[1]+currentCell[2]+currentCell[3]);

        // 4 to 7
        velocity[0] += -currentCell[5]+currentCell[7];
        velocity[1] += currentCell[4]-(currentCell[5]+currentCell[6]+currentCell[7]);
        velocity[2] += -currentCell[4];

        // 8 to 11
        velocity[0] += -currentCell[8]+currentCell[10]-currentCell[11];
        velocity[1] += currentCell[11];

        // 12 to 15
        velocity[0] += currentCell[13]-currentCell[15];
        velocity[1] += currentCell[12]+currentCell[13]-currentCell[14];
        velocity[2] += currentCell[14]+currentCell[15];

        // 16 to 18
        velocity[0] += currentCell[17];
        velocity[1] += currentCell[18];
        velocity[2] += currentCell[16]+currentCell[17]+currentCell[18];

        // Divide by density
        velocity[0] /= (*density);
        velocity[1] /= (*density);
        velocity[2] /= (*density);
    } else {
        velocity[0] = 0.0;
        velocity[1] = 0.0;
        velocity[2] = 0.0;
    }
}

void computeDensityAndVelocity(t_component *c, const t_procData * const procData){

    for (int z = 1; z <= procData->xLength[2] ; z++) {
        for (int y = 1; y <= procData->xLength[1]; y++) {
            for (int x = 1; x <= procData->xLength[0]; x++) {

                double commonVelocity[3];
                double c_velocity[numComp][3];
                double c_density[numComp];
                int fieldIdx = p_computeCellOffsetXYZ(x, y, z, procData->xLength);
                int cellIdx = Q*fieldIdx;

                for(int k = 0; k < numComp; ++k){
                    double *currentCell = &c[k].streamField[cellIdx];

                    // compute density
                    computeNumDensity(currentCell, &c_density[k]);
                    c[k].rho[fieldIdx] = c_density[k];

                    // Compute component velocity
                    computeVelocity(currentCell, &c[k].rho[fieldIdx], c_velocity[k]);
                }
                // Compute common velocity
                computeCommonVelocity(c_density, c_velocity, c, commonVelocity);

                for (int k = 0; k < numComp; k++) {
                    // Compute equilibrium velocity
                    computeEqVelocity(c, commonVelocity, c_density[k], c[k].force[fieldIdx], c[k].velocity[fieldIdx]);
                }
            }
        }
    }
}



/*  Computes the equilibrium distributions for all particle distribution
 *  functions of one cell from density and velocity and stores the results in feq.
 */
void computeFeqCell(const double * const density, const double * const velocity, double *feq){

    // Temporary variables for speed of sound squared and ^4
	// Since it is called that often and having the most work, we made these static
	// to only compute these values once!
    static double const cs_2 = C_S*C_S;
    static double const cs_4_2 = 2*C_S*C_S*C_S*C_S;

    // Temporary variable for velocity
    double const ux = velocity[0];
    double const uy = velocity[1];
    double const uz = velocity[2];

    // Temporary variable for 1-(u.u)/(2*cs²)
    double const u_u = 1-(ux*ux+uy*uy+uz*uz)/(2*cs_2);

    // Temporary variables for density*LATTICEWEIGHTS
    // There are only 3 different LATTTICEWEIGHTS
    double const d1 = (*density)*LATTICEWEIGHTS[0];
    assert(d1 == w3 * (*density));
    double const d2 = (*density)*LATTICEWEIGHTS[2];
    assert(d2 == w2 * (*density));
    double const d3 = (*density)*LATTICEWEIGHTS[9];
    assert(d3 == w1 * (*density));

    // Unroll loop
    // Faster even with -O3
    feq[0]  = d1*(u_u + (-uy-uz)*(1/cs_2 + (-uy-uz)/cs_4_2));
    feq[1]  = d1*(u_u + (-ux-uz)*(1/cs_2 + (-ux-uz)/cs_4_2));
    feq[2]  = d2*(u_u + (-uz)*(1/cs_2 + (-uz)/cs_4_2));
    feq[3]  = d1*(u_u + (ux-uz)*(1/cs_2 + (ux-uz)/cs_4_2));
    feq[4]  = d1*(u_u + (uy-uz)*(1/cs_2 + (uy-uz)/cs_4_2));
    feq[5]  = d1*(u_u + (-ux-uy)*(1/cs_2 + (-ux-uy)/cs_4_2));
    feq[6]  = d2*(u_u + (-uy)*(1/cs_2 + (-uy)/cs_4_2));
    feq[7]  = d1*(u_u + (ux-uy)*(1/cs_2 + (ux-uy)/cs_4_2));
    feq[8]  = d2*(u_u + (-ux)*(1/cs_2 + (-ux)/cs_4_2));
    feq[9]  = d3*(u_u);
    feq[10] = d2*(u_u + (ux)*(1/cs_2 + (ux)/cs_4_2));
    feq[11] = d1*(u_u + (-ux+uy)*(1/cs_2 + (-ux+uy)/cs_4_2));
    feq[12] = d2*(u_u + (uy)*(1/cs_2 + (uy)/cs_4_2));
    feq[13] = d1*(u_u + (ux+uy)*(1/cs_2 + (ux+uy)/cs_4_2));
    feq[14] = d1*(u_u + (-uy+uz)*(1/cs_2 + (-uy+uz)/cs_4_2));
    feq[15] = d1*(u_u + (-ux+uz)*(1/cs_2 + (-ux+uz)/cs_4_2));
    feq[16] = d2*(u_u + (uz)*(1/cs_2 + (uz)/cs_4_2));
    feq[17] = d1*(u_u + (ux+uz)*(1/cs_2 + (ux+uz)/cs_4_2));
    feq[18] = d1*(u_u + (uy+uz)*(1/cs_2 + (uy+uz)/cs_4_2));

    #ifndef NDEBUG
    for (int i = 0; i < Q; i++) {
        assert(feq[i] > 0.0);
    }
    #endif
}

void computeFeq(t_component *c, const t_procData * const procData){

	int cellIdx, fieldIdx;
    for(int k  = 0; k < numComp; ++k){
        for (int z = 1; z <= procData->xLength[2] ; z++) {
            for (int y = 1; y <= procData->xLength[1]; y++) {
                for (int x = 1; x <= procData->xLength[0]; x++) {

                    fieldIdx = p_computeCellOffsetXYZ(x, y, z, procData->xLength);
                    cellIdx = Q*fieldIdx;

                    computeFeqCell(&c[k].rho[fieldIdx], c[k].velocity[fieldIdx], &c[k].feq[cellIdx]);
                }
            }
        }
    }
}


void computeCommonVelocity(const double *const c_density, double c_velocity[2][3], t_component *c, double* commonVel){
    double den = 0.0;
    double momentum[3] = {0.0,0.0,0.0};

    for (int n = 0; n < numComp; ++n) {
       momentum[0]+= c_density[n]*c_velocity[n][0]/c[n].tau;
       momentum[1]+= c_density[n]*c_velocity[n][1]/c[n].tau;
       momentum[2]+= c_density[n]*c_velocity[n][2]/c[n].tau;

       den+= c_density[n]/c[n].tau;
    }
    if (den != 0.0) {
        commonVel[0] = momentum[0]/den;
        commonVel[1] = momentum[1]/den;
        commonVel[2] = momentum[2]/den;
    } else {
        commonVel[0] = 0.0;
        commonVel[1] = 0.0;
        commonVel[2] = 0.0;
    }
}

void computeCellForce(const int currentCellIndex, const int currentCompIndex,
    const t_component *const c, const int * const flagField,
    double const*const G, double forces[3], const t_procData *const procData) {

    // Important: The currentCellIndex is without multiplication with Q
    int xlen2 = procData->xLength[0]+2;
    int xlenylen2 = xlen2*(procData->xLength[1]+2);

    double numDensity;
    forces[0] = 0.0; forces[1] = 0.0; forces[2] = 0.0;

    for (int m = 0; m < numComp; m++) {
        for (int i = 0; i < Q; i++) {

            // Go to the next cell index in the direction of lattice velocities
            int nextCellIndex = currentCellIndex+LATTICEVELOCITIES[i][0]
                                + xlen2*LATTICEVELOCITIES[i][1]
                                + xlenylen2*LATTICEVELOCITIES[i][2]; //index of cell in direction i

            numDensity = c[m].rho[nextCellIndex];

            double G_cur;
            if(LATTICEWEIGHTS[i] == w2){ // 1/18
                G_cur = G[m];
            }else if(LATTICEWEIGHTS[i] == w3){ // 1/36
                G_cur = G[m]/2;
            }else{
                assert(i == 9);
                G_cur = 0.0;
            }

            //Shan&Doolen eq. 4 (PDF page 5)
            forces[0] += G_cur * psiFctPointer[c[m].psiFctCode](numDensity) * LATTICEVELOCITIES[i][0];
            forces[1] += G_cur * psiFctPointer[c[m].psiFctCode](numDensity) * LATTICEVELOCITIES[i][1];
            forces[2] += G_cur * psiFctPointer[c[m].psiFctCode](numDensity) * LATTICEVELOCITIES[i][2];
         }
    }

    numDensity = c[currentCompIndex].rho[currentCellIndex];

    forces[0] *= -psiFctPointer[c[currentCompIndex].psiFctCode](numDensity);
    forces[1] *= -psiFctPointer[c[currentCompIndex].psiFctCode](numDensity);
    forces[2] *= -psiFctPointer[c[currentCompIndex].psiFctCode](numDensity);
}

void computeForce(t_component *c, const t_procData * const procData, int *flagField, double G[numComp][numComp]){

    for(int k = 0; k < numComp; ++k){
        for (int z = 1; z <= procData->xLength[2] ; z++) {
            for (int y = 1; y <= procData->xLength[1]; y++) {
                for (int x = 1; x <= procData->xLength[0]; x++) {
                    int fieldIdx = p_computeCellOffsetXYZ(x, y, z, procData->xLength);
                    computeCellForce(fieldIdx, k, c, flagField, G[k],  c[k].force[fieldIdx], procData);
                }
            }
        }
    }

}

void computeEqVelocity(t_component const*const c, double const*const commonVelocity, const double compDensity, double const*const compForce, double compEqVelocity[3]) {
    if (compDensity != 0.0) {
        compEqVelocity[0] = commonVelocity[0] + (c->tau/compDensity)*compForce[0];
        compEqVelocity[1] = commonVelocity[1] + (c->tau/compDensity)*compForce[1];
        compEqVelocity[2] = commonVelocity[2] + (c->tau/compDensity)*compForce[2];
    } else {
        compEqVelocity[0] = 0.0;
        compEqVelocity[1] = 0.0;
        compEqVelocity[2] = 0.0;
    }
}
