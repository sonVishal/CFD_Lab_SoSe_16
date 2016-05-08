#include <stdio.h>
#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

/* TODO: (DL)
 * NOTE IF SOMEONE WANTS TO CONTINUE HERE:
 *
 * I'm not sure if this is the 'best' way how I approached this in the first go (got more complex than I first thought).
 *
 * The basic idea is, to not go through ALL cells but only the boundary cells to reduce read accesses.
 *
 * So for each fixed (x,y,z) there is one function at the moment -- the 'direction' indicates in which direction
 * the domain lies (e.g. x is fixed, +1 -> in positive x direction there are in-domain cells).
 *
 * At the moment there is a lot of redundance between these pWall functions - they have only different index
 * computations. Maybe later on it is worth to handle callbacks and only have 1 function.
 *
 * At the moment dont trust the indices etc. too much - until now the main importance was the structure!
 *
 */

/* TODO: (DL)
 * Needs to be static when using Q and C_S - alternative would be to handle the values
 * (and make function not static
 */

/* TODO (DL) Think about: join all these pWall functions and provide callback functions that determine how
 * to compute the index!
 */

static inline void setBounceBack(double *collideField, const double * const wallVelocity,
		const int type, const int i, const int cell_xyzoffset, const int n_xyzoffset, const int *c) {

	if(type == 1){
		collideField[cell_xyzoffset + i] = collideField[n_xyzoffset + (Q-i-1)];
	}else if(type == 2){
		//c_2 = ||c_i||^2 -- squared to need not take the square root
		int c_2 = c[0] + c[1] + c[2];

		if(c_2 == 0){
			printf("WARNING: c_2 is zero, which should not happen here!!! \n");
		}

		double w_i = c_2 == 1 ? w2 : w3;

		double density = 0;
		double dot_uwall_c = 0;
		computeDensity(&collideField[n_xyzoffset], &density);
		pDotProduct2(wallVelocity, &c[0], &dot_uwall_c);

		collideField[cell_xyzoffset + i] = collideField[n_xyzoffset + (Q-i-1)] +
				2 * w_i * dot_uwall_c / (C_S*C_S);

	}else{
		printf("WARNING: a FLUID cell appeared - this should not happen!!\n");
	}
}

void pxWalls(double *collideField, const int * const flagField, const double * const wallVelocity, const int x,
		const int xlength, const int direction){

	int xlength_2 = xlength*xlength;

	for(int z = 0; z <= xlength; ++z){
		int zoffset = z*xlength_2;
		for(int y = 0; y <= xlength; ++y){
			int xyz_offset = zoffset + y*xlength + x;
			int cell_xyzoffset = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				if(LATTICEVELOCITIES[i][0] == direction){
					int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
					int n_xyzoffset = Q*( (z+c[0])*xlength_2 + (y+c[1])*xlength + x+c[0] );
					setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, cell_xyzoffset, n_xyzoffset, c);
				}
			}
		}
	}
}

void pyWalls(double *collideField, const int * const flagField, const double * const wallVelocity, const int y,
		const int xlength, const int direction){
	int xlength_2 = xlength*xlength;
	int yoffset = y*xlength;

	for(int z = 0; z < xlength; ++z){
		int zoffset = z*xlength_2;
		for(int x = 0; x < xlength; ++x){

			int xyz_offset = zoffset + yoffset + x;
			int cell_xyzoffset = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				if(LATTICEVELOCITIES[i][1] == direction){
					int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
					int n_xyzoffset = Q*( (z+c[0])*xlength_2 + (y+c[1])*xlength + x+c[0] );
					setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, cell_xyzoffset, n_xyzoffset, c);
				}
			}
		}
	}
}

void pzWalls(double *collideField, const int * const flagField, const double * const wallVelocity, const int z,
		const int xlength, const int direction){
	int xlength_2 = xlength * xlength;
	int zoffset = z*xlength*xlength;

	for(int y = 0; y < xlength; ++y){
		int yzoffset = y*xlength + zoffset;
		for(int x = 0; x < xlength; ++x){

			int xyz_offset = yzoffset + x;
			int cell_xyzoffset = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				if(LATTICEVELOCITIES[i][2] == direction){
					int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
					int n_xyzoffset = Q*( (z+c[0])*xlength_2 + (y+c[1])*xlength + x+c[0] );
					setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, cell_xyzoffset, n_xyzoffset, c);
				}
			}
		}
	}
}

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
	//Loop over boundary cells only

	/* TODO: Try to find optimizing ways - not cache friendly walls: */
	pxWalls(collideField, flagField, wallVelocity, 0,       xlength,  1);
	pxWalls(collideField, flagField, wallVelocity, xlength, xlength, -1);

	pyWalls(collideField, flagField, wallVelocity, 0,       xlength,  1);
	pyWalls(collideField, flagField, wallVelocity, xlength, xlength, -1);

	/* most cache friendly walls */
	pzWalls(collideField, flagField, wallVelocity, 0,       xlength,  1);
	pzWalls(collideField, flagField, wallVelocity, xlength, xlength, -1);
}
