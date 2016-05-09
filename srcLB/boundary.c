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
 * Also I wanted to implement it cache efficient. That is, the different boundary walls are
 * treated separately (to prevent too many scattered accesses).
 *
 * This way the z-fixed fall and the y-fixed boundary walls profit from cache. The x-fixed wall is bad.
 *
 * So for each fixed (x,y,z) there is one function at the moment -- the 'direction' indicates in which direction
 * the domain lies (e.g. x is fixed, +1 -> in positive x direction there are in-domain cells).
 *
 * At the moment there is a lot of redundance between these pWall functions - they have only different index
 * computations. Maybe later on it is worth to handle callbacks and only have 1 function.
 *
 * At the moment dont trust the indices etc. too much - until now the main importance was the structure!
 */

/* TODO: (DL)
 * Needs to be static when using Q and C_S - alternative would be to handle the values
 * (and make function not static
 */

/* TODO (DL) Think about: join all these pWall functions and provide callback functions that determine how
 * to compute the index!
 */

static inline void setBounceBack(double *collideField, const double * const wallVelocity,
		const int type, const int i, const int current_cell_index, const int n_cell_index, const int *c) {

	if(type == 1){
		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)];
	}else if(type == 2){

		double density = 0;
		double dot_uwall_c = 0;
		computeDensity(&collideField[n_cell_index], &density);
		pDotProduct2(wallVelocity, &c[0], &dot_uwall_c);

		double weight = LATTICEWEIGHTS[i];

		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)] +
				2 * weight * dot_uwall_c / (C_S*C_S);

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
			int current_cell_index = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				if(LATTICEVELOCITIES[i][0] == direction){
					int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
					int n_cell_index = Q*( (z+c[0])*xlength_2 + (y+c[1])*xlength + x+c[0] );
					setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, current_cell_index, n_cell_index, c);
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
			int current_cell_index = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				if(LATTICEVELOCITIES[i][1] == direction){
					int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
					int n_cell_index = Q*( (z+c[0])*xlength_2 + (y+c[1])*xlength + x+c[0] );
					setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, current_cell_index, n_cell_index, c);
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
			int current_cell_index = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				if(LATTICEVELOCITIES[i][2] == direction){
					int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
					int n_cell_index = Q*( (z+c[0])*xlength_2 + (y+c[1])*xlength + x+c[0] );
					setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, current_cell_index, n_cell_index, c);
				}
			}
		}
	}
}

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
	//Loop over boundary cells only

	/* TODO: Try to find optimizing ways - not cache friendly walls: */
	pxWalls(collideField, flagField, 0, 0,       xlength,  1);
	pxWalls(collideField, flagField, 0, xlength + 1, xlength, -1);

	pyWalls(collideField, flagField, 0, 0,       xlength,  1);
	pyWalls(collideField, flagField, 0, xlength + 1, xlength, -1);

	/* most cache friendly walls */
	pzWalls(collideField, flagField, 0, 0,       xlength,  1);
	pzWalls(collideField, flagField, wallVelocity, xlength + 1, xlength, -1);
}
