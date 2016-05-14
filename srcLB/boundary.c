#include <stdio.h>
#include "boundary.h"
#include "helper.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

/*
 * TODO: (DL) treat overlapping cells (double computation of corner values)
 */

/* Enum that gives a distinction of the walls with respective to the fixed value of (x,y,z) */
enum WALLS{E_XFIXED, E_YFIXED, E_ZFIXED};

/* Helper function that carries out the moving wall or bounce back condition, depending on flag 'type'.
 * The value in flagField corresponding to the current cell (current_cell_index) has to indicate a boundary
 * cell (flag 1 or 2), the value in flagField corresponding to neighboring cell (n_cell_index) has to
 * be a cell in the domain (flag 0).
 */
void p_setBounceBack(double *collideField, const double * const wallVelocity,
		const int type, const int i, const int current_cell_index, const int n_cell_index, const int *c) {

	if(type == 1){ 			//bounce back
		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)];
	}else if(type == 2){ 	//moving wall

		double density = 0;
		computeDensity(&collideField[n_cell_index], &density);

		double dot_uwall_c = wallVelocity[0]*c[0]+wallVelocity[1]*c[1]+wallVelocity[2]*c[2];
		double weight = LATTICEWEIGHTS[i];

		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)] +
				2 * weight * dot_uwall_c / (C_S*C_S);

	}else{ 		//Wrong input parameter
		ERROR("A FLUID cell appeared to set boundaries. This should not happen!!\n");
	}
}

/* Helper function that computes the offset of the current cell. 'Inner' corresponds to the first value of (x,y,z)
 * that is not fixed, 'outer' to the second value. By this a better cache efficiency is obtained.
 * WallIdx has to be a valid index from the enum WALLS.
 */
int p_computeCellOffset(const int outer, const int inner, const int fixedValue, const int xlength, const int wallIdx){
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
			ERROR("Invalid wall index occured. This should not happen !!!");
			return -1;
	}
}

/* Helper function similar to p_computeCellOffset, except that it computes the offset value
 * of a neighboring cell in the direction of 'c_vec'.
 * There are no checks that guarantee that the computed offset value is valid.
 * For example, when the current cell is on the boundary and 'c_vec' points away from the inner value.
 */
int p_computeNeighborCellOffset(int outer, int inner, int fixedValue,
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
			ERROR("Invalid wall index occured. This should not happen !!!");
			return -1;
	}
}

/* Helper function that treats a single boundary wall. The wall itself is described with wallIdx
 * (from enum WALLS) and the 'fixedValue' (generally 0 or xlength+1).
 */

void p_treatSingleWall(double *collideField, const int * const flagField,
		const int fixedValue, const double * const wallVelocity, const int xlength, const int wallIdx){

	//needed to check whether it is a potential in-domain cell
	int maxValidIndex = Q*(xlength+2)*(xlength+2)*(xlength+2);

	//k - corresponds to the 'outer' value when computing the offset
	for(int k = 0; k <= xlength+1; ++k){
		//j - corresponds to the 'inner' value
		for(int j = 0; j <= xlength+1; ++j){

			//flagField[xyz_offset] should only be a boundary cell (checked in p_setBounceBack)
			int xyz_offset = p_computeCellOffset(k, j, fixedValue, xlength, wallIdx);
			int current_cell_index = Q*xyz_offset;

			for(int i = 0; i < Q; ++i){
				int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
				int n_xyzoffset = p_computeNeighborCellOffset(k, j, fixedValue, c, xlength, wallIdx);
				int n_cell_index = Q*n_xyzoffset;

				//check valid index: in case the direction of vector 'c' points to a non-existing cell
				//check that the neighboring cell is a FLUID cell in the domain (and not another boundary cell)
				if(n_cell_index >= 0 && n_cell_index < maxValidIndex &&
						flagField[n_xyzoffset] == 0
				){
					p_setBounceBack(collideField, wallVelocity, flagField[xyz_offset], i, current_cell_index, n_cell_index, c);
				}
			}
		}
	}
}

/* Helper function that does both boundary walls for the fixed value provided in the enum WALLS.
 * The fixed value is set to 0 and xlength+1.
 */
void p_setWallBoundaries(double *collideField, const int * const flagField,
		const double * const wallVelocity, const int xlength, const int wallIdx){

	//Fixed value of 0
	p_treatSingleWall(collideField, flagField, 0,         wallVelocity, xlength, wallIdx);

	//Fixed value of xlength+1
	p_treatSingleWall(collideField, flagField, xlength+1, wallVelocity, xlength, wallIdx);
}

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
	p_setWallBoundaries(collideField, flagField, wallVelocity, xlength, E_XFIXED);
	p_setWallBoundaries(collideField, flagField, wallVelocity, xlength, E_YFIXED);
	p_setWallBoundaries(collideField, flagField, wallVelocity, xlength, E_ZFIXED);
}
