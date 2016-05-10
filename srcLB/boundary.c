#include <stdio.h>
#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

/*
 * TODO: (DL) treat overlapping cells (double computation of corner values)
 */

enum WALLS{
	E_XFIXED,
	E_YFIXED,
	E_ZFIXED
};

void setBounceBack(double *collideField, const double * const wallVelocity,
		const int type, const int i, const int current_cell_index, const int n_cell_index, const int *c) {

	if(type == 1){
		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)];
	}else if(type == 2){

		double density = 0;
		computeDensity(&collideField[n_cell_index], &density);

		double dot_uwall_c = wallVelocity[0]*c[0]+wallVelocity[1]*c[1]+wallVelocity[2]*c[2];
		double weight = LATTICEWEIGHTS[i];

		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)] +
				2 * weight * dot_uwall_c / (C_S*C_S);

	}else{
		printf("WARNING: a FLUID cell appeared - this should not happen!!\n");
	}
}

int computeCellOffset(const int outer, const int inner, const int fixedValue, const int xlength, const int wallIdx){
	switch (wallIdx) {
		case E_XFIXED:
			//outer = Z, inner = Y
			return (xlength+2) * (outer*(xlength+2) + inner) + fixedValue;
		case E_YFIXED:
			//outer = Z, inner = X
			return (xlength+2) * (outer*(xlength+2) + fixedValue) + inner;
		case E_ZFIXED:
			//outer = Y, inner = X
			return (xlength+2) * (fixedValue*(xlength+2) + outer) + inner;
		default:
			printf("ERROR: INVALID WALL INDEX OCCURED. THIS SHOULD NOT HAPPEN!!!");
			return -1;
	}
}

int computeNeighborCellOffset(int outer, int inner, int fixedValue,
		const int * const c_vec, const int xlength, const int wallIdx){
	switch (wallIdx) {
		case E_XFIXED:
			outer += c_vec[2];      //= Z
			inner += c_vec[1];      //= Y
			fixedValue += c_vec[0]; //= X
			return (xlength+2) * (outer*(xlength+2) + inner) + fixedValue;
		case E_YFIXED:
			outer += c_vec[2];      //= Z
			inner += c_vec[0];      //= X
			fixedValue += c_vec[1]; //= Y
			return (xlength+2) * (outer*(xlength+2) + fixedValue) + inner;
		case E_ZFIXED:
			outer += c_vec[1];      //= Y
			inner += c_vec[0];      //= X
			fixedValue += c_vec[2]; //= Z
			return (xlength+2) * (fixedValue*(xlength+2) + outer) + inner;

		default:
			printf("ERROR: INVALID WALL INDEX OCCURED. THIS SHOULD NOT HAPPEN!!!");
			return -1;
	}
}

void treatSingleWall(double *collideField, const int * const flagField,
		const int fixedValue, const double * const wallVelocity, const int xlength, const int wallIdx){

	//needed to check whether it is a potential in-domain cell
	int maxValidIndex = Q*(xlength+2)*(xlength+2)*(xlength+2);

	for(int k = 0; k <= xlength+1; ++k){
		for(int j = 0; j <= xlength+1; ++j){
			int xyz_offset = computeCellOffset(k, j, fixedValue, xlength, wallIdx);
			int current_cell_index = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
				int n_xyzoffset = computeNeighborCellOffset(k, j, fixedValue, c, xlength, wallIdx);
				int n_cell_index = Q*n_xyzoffset;

				if(n_cell_index >= 0 && n_cell_index < maxValidIndex && //check valid index
						flagField[n_xyzoffset] == 0 // check if neighbor is FLUID field (and not another boundary cell)
				){
					setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, current_cell_index, n_cell_index, c);
				}
			}
		}
	}
}

void setWallBoundaries(double *collideField, const int * const flagField,
		const double * const wallVelocity, const int xlength, const int wallIdx){

	//Fixed value of 0
	treatSingleWall(collideField, flagField, 0,         wallVelocity, xlength, wallIdx);

	//Fixed value of xlength+1
	treatSingleWall(collideField, flagField, xlength+1, wallVelocity, xlength, wallIdx);
}


void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
	//Loop over boundary cells only
	setWallBoundaries(collideField, flagField, wallVelocity, xlength, E_XFIXED);
	setWallBoundaries(collideField, flagField, wallVelocity, xlength, E_YFIXED);
	setWallBoundaries(collideField, flagField, wallVelocity, xlength, E_ZFIXED);
}
