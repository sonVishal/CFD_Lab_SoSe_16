#include "computeCellValues.h"
#include "LBDefinitions.h"

//TODO: (TKS) Discuss convention that component specific functions hav c_ prefix


/** computes the number density and density from the particle distribution 
 * functions stored at  currentCell. currentCell thus denotes the address 
 * of the first particle  distribution function of the respective cell.
 *  The result is stored in density.
 */


void c_computeNumDensity(const double *const currentCell, double *c_numDensity){
    // Number density is the sum of the distributions in the current lattice
    int i;
    *c_numDensity = 0.0;
    for (i = 0; i < Q; i++) {
		 *c_numDensity += currentCell[i];
    }
}


//TODO: (TKS) Not changed yet
/** computes the velocity within currentCell and stores the result in velocity */
void c_computeVelocity(const double * const currentCell, const double * const density, double *velocity, const double * const mass){

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

    // Divide by density and multiply by mass
    velocity[0] *= *mass/(*density);
    velocity[1] *= *mass/(*density);
    velocity[2] *= *mass/(*density);
}

/*computes the total velociy as if there were no interacting forces*/
void computeVelocityNI(const int* numComp, const double *const c_density, const double * const c_velocity, const double *const c_tau, double* velocityNI){

    double den = 0;
    double nom[3] = {0,0,0};

    for (int n = 0; n < *numComp; ++n) {
       nom[0]+= c_density[n]*c_velocity[0]/c_tau[n];
       nom[1]+= c_density[n]*c_velocity[1]/c_tau[n];
       nom[2]+= c_density[n]*c_velocity[2]/c_tau[n];

       den+= c_density[n]/c_tau[n];
    }

    velocityNI[0] = nom[0]/den;
    velocityNI[1] = nom[1]/den;
    velocityNI[2] = nom[2]/den;

}

/*computes interacting forces between species*/
//
//TODO: (TKS) Add const where appropriate
//TODO: (TKS) Tidy up signature
void computeForces(int currentCellIndex, t_component *c, const int *numComp, double *c_numDensity, double **G, int * xlength, double *forces[3]){
    //TODO: (TKS) Take into consideration only neighbours.

    int xlen2 = xlength[0]+2;
	int ylen2 = xlength[1]+2;
    int xylen = xlen2*ylen2;

    double numDensity;
    
    for (int n = 0; n < *numComp; ++n){
        forces[n] = 0;

        for (int m = 0; m < *numComp; ++m) {
            
            //TODO: (TKS) Could unroll loop.
            for (int i = 0; i < Q; i++) {
                int nextCellIndex = 
                currentCellIndex-LATTICEVELOCITIES[i][0]-
                xlen2*LATTICEVELOCITIES[i][1]-
                xylen*LATTICEVELOCITIES[i][2]; //index of cell in direction i

                int nextIndex = Q*nextCellIndex; //index of number density in direction i
                numDensity = c[m].collideField[nextIndex]; //number density in direction i


                forces[n][0] += G[n][m] * c[m].psi(numDensity) * LATTICEVELOCITIES[i][0];
                forces[n][1] += G[n][m] * c[m].psi(numDensity) * LATTICEVELOCITIES[i][1];
                forces[n][2] += G[n][m] * c[m].psi(numDensity) * LATTICEVELOCITIES[i][2];
             }

        }

        numDensity = c[n].collideField[currentCellIndex];
        forces[n][0] *= c[n].psi(numDensity);
        forces[n][1] *= c[n].psi(numDensity);
        forces[n][2] *= c[n].psi(numDensity);
    }
}

/** computes the equilibrium distributions for all particle distribution
 *  functions of one cell from density and velocity and stores the results in feq.
 */
void computeFeq(const double density, const double * const velocity, double *feq){

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
    /*double const u_u = 1-(ux*ux+uy*uy+uz*uz)/(2*cs_2);*/

    // Temporary variables for density*LATTICEWEIGHTS
    // There are only 3 different LATTTICEWEIGHTS
    /*double const d1 = density*LATTICEWEIGHTS[0];*/
    /*double const d2 = density*LATTICEWEIGHTS[2];*/
    /*double const d3 = density*LATTICEWEIGHTS[9];*/

    // Semantics for the unrolled loop
    // int i;
    //
    // for (i = 0; i < Q; i++) {
	// 	double dotProd1 = ux*ux + uy*uy + uz*uz;
	// 	double dotProd2 = ux*LATTICEVELOCITIES[i][0]
    //                     + uy*LATTICEVELOCITIES[i][1]
    //                     + uz*LATTICEVELOCITIES[i][2];
	// 	feq[i] = LATTICEWEIGHTS[i]*(*density)*(1 + dotProd2/cs_2 + dotProd2*dotProd2/cs_4_2 - dotProd1/(2*cs_2));
	// }

    // Unroll loop
    // Faster even with -O3
    //feq[0]  = d1*(u_u + (-uy-uz)*(1/cs_2 + (-uy-uz)/cs_4_2));
    //feq[1]  = d1*(u_u + (-ux-uz)*(1/cs_2 + (-ux-uz)/cs_4_2));
    //feq[2]  = d2*(u_u + (-uz)*(1/cs_2 + (-uz)/cs_4_2));
    //feq[3]  = d1*(u_u + (ux-uz)*(1/cs_2 + (ux-uz)/cs_4_2));
    //feq[4]  = d1*(u_u + (uy-uz)*(1/cs_2 + (uy-uz)/cs_4_2));
    //feq[5]  = d1*(u_u + (-ux-uy)*(1/cs_2 + (-ux-uy)/cs_4_2));
    //feq[6]  = d2*(u_u + (-uy)*(1/cs_2 + (-uy)/cs_4_2));
    //feq[7]  = d1*(u_u + (ux-uy)*(1/cs_2 + (ux-uy)/cs_4_2));
    //feq[8]  = d2*(u_u + (-ux)*(1/cs_2 + (-ux)/cs_4_2));
    //feq[9]  = d3*(u_u);
    //feq[10] = d2*(u_u + (ux)*(1/cs_2 + (ux)/cs_4_2));
    //feq[11] = d1*(u_u + (-ux+uy)*(1/cs_2 + (-ux+uy)/cs_4_2));
    //feq[12] = d2*(u_u + (uy)*(1/cs_2 + (uy)/cs_4_2));
    //feq[13] = d1*(u_u + (ux+uy)*(1/cs_2 + (ux+uy)/cs_4_2));
    //feq[14] = d1*(u_u + (-uy+uz)*(1/cs_2 + (-uy+uz)/cs_4_2));
    //feq[15] = d1*(u_u + (-ux+uz)*(1/cs_2 + (-ux+uz)/cs_4_2));
    //feq[16] = d2*(u_u + (uz)*(1/cs_2 + (uz)/cs_4_2));
    //feq[17] = d1*(u_u + (ux+uz)*(1/cs_2 + (ux+uz)/cs_4_2));
    //feq[18] = d1*(u_u + (uy+uz)*(1/cs_2 + (uy+uz)/cs_4_2));
    
    //TODO: (TKS) Unroll the loop. Much easier to debug this.
     int i;
    
     for (i = 0; i < Q; i++) {
         double dotProd1 = ux*ux + uy*uy + uz*uz;
         double dotProd2 = ux*LATTICEVELOCITIES[i][0]
                         + uy*LATTICEVELOCITIES[i][1]
                         + uz*LATTICEVELOCITIES[i][2];
         feq[i] = LATTICEWEIGHTS[i]*(density)*(1 + dotProd2/cs_2 + dotProd2*dotProd2/cs_4_2 - dotProd1/(2*cs_2));
     }

}
