#include "computeCellValues.h"

/** computes the density from the particle distribution functions stored at
 *  currentCell. currentCell thus denotes the address of the first particle
 *  distribution function of the respective cell.
 *  The result is stored in density.
 */
void computeDensity(const double *const currentCell, double *density){
    /* Semantics of unrolled loop:
     *
	 *	for (i = 0; i < Q; i++) {
	 *		 *density += currentCell[i];
	 *	}
     */

    (*density) = currentCell[0]+currentCell[1]+currentCell[2]+currentCell[3]+
        currentCell[4]+currentCell[5]+currentCell[6]+currentCell[7]+
        currentCell[8]+currentCell[9]+currentCell[10]+currentCell[11]+
        currentCell[12]+currentCell[13]+currentCell[14]+currentCell[15]+
        currentCell[16]+currentCell[17]+currentCell[18];
}

/** computes the velocity within currentCell and stores the result in velocity */
void computeVelocity(const double * const currentCell, const double * const density, double *velocity){

	/* Semantics of unrolled loop:
	 *
	 * for (i = 0; i < Q; i++) {
	 * 	   velocity[0] += currentCell[i]*LATTICEVELOCITIES[i][0];
	 * 	   velocity[1] += currentCell[i]*LATTICEVELOCITIES[i][1];
	 * 	   velocity[2] += currentCell[i]*LATTICEVELOCITIES[i][2];
	 * }
	 */

	// Unroll the loop for cache efficient access
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
}

/** computes the equilibrium distributions for all particle distribution
 *  functions of one cell from density and velocity and stores the results in feq.
 */
void computeFeq(const double * const density, const double * const velocity, double *feq){

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
    double const d2 = (*density)*LATTICEWEIGHTS[2];
    double const d3 = (*density)*LATTICEWEIGHTS[9];


    /* Semantics of unrolled loop:
     * for (i = 0; i < Q; i++) {
	 *		dot_product(velocity, velocity);
	 *		dot_product(velocity, LATTICEVELOCITIES[i]);
	 *		tmp = dotProd2/C_S_2;
	 *		feq[i] = LATTICEWEIGHTS[i]*(*density)*(1 + tmp + tmp*tmp/2 - dotProd1/(2*C_S_2));
	 *	}
	 */

    // Unroll loop
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

}
