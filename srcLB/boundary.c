#include <stdio.h>
#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

/* TODO: (DL)
 * Needs to be static when using Q and C_S - alternative would be to handle the values
 * (and make function not static
 */

/* TODO (DL) Think about: somehow join these pWall functions
 */

/*
 * TODO: (DL) treat overlapping cells (double computation of corner values)
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

	int xlength_2 = (xlength+2)*(xlength+2);

	for(int z = 0; z <= xlength+1; ++z){
		int zoffset = z*xlength_2;
		for(int y = 0; y <= xlength+1; ++y){
			int xyz_offset = zoffset + y*(xlength+2) + x;
			int current_cell_index = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
				int n_xyzoffset = (z+c[2])*xlength_2 + (y+c[1])*(xlength+2) + x+c[0];
				int n_cell_index = Q*n_xyzoffset;
				if(n_cell_index >= 0 && n_cell_index <=Q*(xlength+2)*(xlength+2)*(xlength+2) && //check valid index
						flagField[n_xyzoffset] == 0 // check if neighbor is FLUID field (and not another boundary cell)
				){
					setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, current_cell_index, n_cell_index, c);
				}
			}
		}
	}
}

void pyWalls(double *collideField, const int * const flagField, const double * const wallVelocity, const int y,
		const int xlength, const int direction){
	int xlength_2 = (xlength+2)*(xlength+2);
	int yoffset = y*(xlength+2);

	for(int z = 0; z <= xlength+1; ++z){
		int zoffset = z*xlength_2;
		for(int x = 0; x <= xlength+1; ++x){

			int xyz_offset = zoffset + yoffset + x;
			int current_cell_index = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
				int n_xyzoffset = (z+c[2])*xlength_2 + (y+c[1])*(xlength+2) + x+c[0];
				int n_cell_index = Q*n_xyzoffset;
				if(n_cell_index >= 0 && n_cell_index <=Q*(xlength+2)*(xlength+2)*(xlength+2) && //check valid index
						flagField[n_xyzoffset] == 0 // check if neighbor is FLUID field (and not another boundary cell)
				){
					setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, current_cell_index, n_cell_index, c);
				}
			}
		}
	}
}

void pzWalls(double *collideField, const int * const flagField, const double * const wallVelocity, const int z,
		const int xlength, const int direction){
	int xlength_2 = (xlength+2) * (xlength+2);
	int zoffset = z*xlength_2;

	for(int y = 0; y <= xlength+1; ++y){
		int yzoffset = y*(xlength+2) + zoffset;
		for(int x = 0; x <= xlength+1; ++x){

			int xyz_offset = yzoffset + x;
			int current_cell_index = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
				int n_xyzoffset = (z+c[2])*xlength_2 + (y+c[1])*(xlength+2) + x+c[0];
				int n_cell_index = Q*n_xyzoffset;
				if(n_cell_index >= 0 && n_cell_index <=Q*(xlength+2)*(xlength+2)*(xlength+2) && //check valid index
						flagField[n_xyzoffset] == 0 // check if neighbor is FLUID field (and not another boundary cell)
				){
					// printf("FLAG ID: %i \n", flagField[xyz_offset]); //TODO: (DL) delete... just for testing
					setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, current_cell_index, n_cell_index, c);
				}
			}
		}
	}
}

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
	//Loop over boundary cells only

	/* most cache friendly walls */
	pzWalls(collideField, flagField, wallVelocity, 0,         xlength,  1);
	pzWalls(collideField, flagField, wallVelocity, xlength+1, xlength, -1);

	/* TODO: Try to find optimizing ways - not cache friendly walls: */
	pxWalls(collideField, flagField, wallVelocity, 0,         xlength,  1);
	pxWalls(collideField, flagField, wallVelocity, xlength+1, xlength, -1);

	pyWalls(collideField, flagField, wallVelocity, 0,         xlength,  1);
	pyWalls(collideField, flagField, wallVelocity, xlength+1, xlength, -1);

}
