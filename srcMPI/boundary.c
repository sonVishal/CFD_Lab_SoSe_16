#include <stdio.h>
#include "boundary.h"
#include "helper.h"
#include "debug.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <mpi/mpi.h>

/* Helper function that carries out the moving wall or bounce back condition, depending on flag 'type'.
 * The value in flagField corresponding to the current cell (current_cell_index) has to indicate a boundary
 * cell (flag 1 or 2), the value in flagField corresponding to neighboring cell (n_cell_index) has to
 * be a cell in the domain (= FLUID cell of flag 0).
 */
void p_setBounceBack(double *collideField, const double * const wallVelocity,
		const int type, const int i, const int current_cell_index, const int n_cell_index, const int *c) {

	if(type == 1){ 			//bounce back
		/* equation 17 */
		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)];
	}else{ // type == 2 (MOVING_WALL)

		double density = 0;
		computeDensity(&collideField[n_cell_index], &density);

		double dot_uwall_c = wallVelocity[0]*c[0]+wallVelocity[1]*c[1]+wallVelocity[2]*c[2];
		double weight = LATTICEWEIGHTS[i];

		/* equation 19 */
		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)] +
				2 * weight * density * dot_uwall_c / (C_S*C_S);

	}
}

/* Helper function that computes the offset of the current cell. 'Inner' corresponds to the first value of (x,y,z)
 * that is not fixed; 'outer' to the second value. By this ordering a better cache efficiency is obtained.
 * WallIdx has to be a valid index from the enum WALLS.
 */
int p_computeCellOffset(const int outer, const int inner, const int fixedValue, int const * const xlength, const int wallIdx){

	//Computes index: z * (xlen*ylen) + y * (xlen) + x

	//wallIdx has valid integer values from 0 to 5
	switch (wallIdx/2) { //integer division to get the type of face (see enum in LBDefinitions.h)
		case 0: // LEFT, RIGHT -> Y fixed
			//outer = Z, inner = X
			return (xlength[0]+2) * (outer * (xlength[1]+2) + fixedValue) + inner;

		case 1: // TOP, BOTTOM -> Z fixed
			//outer = Y, inner = X
			return (xlength[0]+2) * (fixedValue * (xlength[1]+2) + outer) + inner;

		case 2: // FRONT, BACK -> X fixed
			//outer = Z, inner = Y
			return (xlength[0]+2) * (outer * (xlength[1]+2) + inner) + fixedValue;

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
		const int * const c_vec, int const * const xlength, const int wallIdx){

	switch (wallIdx/2) {
		case 0: // LEFT, RIGHT -> Y fixed
			outer 		+= c_vec[2];	//= Z
			inner 		+= c_vec[0];	//= X
			fixedValue 	+= c_vec[1];	//= Y
			return (xlength[0]+2) * (outer * (xlength[1]+2) + fixedValue) + inner;

		case 1: // TOP, BOTTOM -> Z fixed
			outer 		+= c_vec[1];	//= Y
			inner 		+= c_vec[0];	//= X
			fixedValue 	+= c_vec[2];	//= Z
			return (xlength[0]+2) * (fixedValue * (xlength[1]+2) + outer) + inner;

		case 2: // FRONT, BACK -> X fixed
			outer 		+= c_vec[2];	//= Z
			inner 		+= c_vec[1];	//= Y
			fixedValue 	+= c_vec[0];	//= X
			return (xlength[0]+2) * (outer * (xlength[1]+2) + inner) + fixedValue;

		default:
			ERROR("Invalid wall index occured. This should not happen !!!");
			return -1;
	}
}


void p_setIterationParameters(int *endOuter, int *endInner, int *fixedValue, const t_procData procData, const int wallIdx){

	switch(wallIdx){

	//---------------------------------------------
	//outer = Z, inner = X, Y fixed
	case LEFT:
		*endOuter = procData.xLength[2]+1;
		*endInner = procData.xLength[0]+1;
		*fixedValue = 0;
		break;
	case RIGHT:
		*endOuter = procData.xLength[2]+1;
		*endInner = procData.xLength[0]+1;
		*fixedValue = procData.xLength[1]+1;
		break;
	//---------------------------------------------
	//outer = Y, inner = X, Z fixed
	case TOP:
		*endOuter = procData.xLength[1]+1;
		*endInner = procData.xLength[0]+1;
		*fixedValue = procData.xLength[2]+1;
		break;
	case BOTTOM:
		*endOuter = procData.xLength[1]+1;
		*endInner = procData.xLength[0]+1;
		*fixedValue = 0;
		break;

	//---------------------------------------------
	//outer = Z, inner = Y, X fixed
	case FRONT:
		*endOuter = procData.xLength[2]+1;
		*endInner = procData.xLength[1]+1;
		*fixedValue = procData.xLength[0]+1;
		break;
	case BACK:
		*endOuter = procData.xLength[2]+1;
		*endInner = procData.xLength[1]+1;
		*fixedValue = 0;
		break;

	default:
		ERROR("Invalid wallIdx occurred. This should never happen!");
	}

}


/* Helper function that treats a single boundary wall. The wall itself is described with wallIdx
 * (from enum WALLS) and the 'fixedValue' (generally 0 or xlength+1).
 */

void p_treatSingleWall(double *collideField, const int * const flagField, const t_procData procData, const int wallIdx){

	int endOuter = -1, endInner = -1, fixedValue = -1;
	p_setIterationParameters(&endOuter, &endInner, &fixedValue, procData, wallIdx);
	assert(endOuter != -1 && endInner != -1 && fixedValue != -1);

	int xlength[3] = {procData.xLength[0], procData.xLength[1], procData.xLength[2]};

	//variable needed at check whether it is an in-domain (FLUID) cell
	int maxValidIndex = Q*(xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2);

	//k - corresponds to the 'outer' value when computing the offset
	for(int k = 0; k <= endOuter; ++k){
		//j - corresponds to the 'inner' value
		for(int j = 0; j <= endInner; ++j){

			//flagField[xyz_offset] should only be a boundary cell (checked in p_setBounceBack)
			int xyz_offset = p_computeCellOffset(k, j, fixedValue, xlength, wallIdx);
			int current_cell_index = Q*xyz_offset;
			int boundaryType = flagField[xyz_offset];

			assert(boundaryType == NO_SLIP || boundaryType == PARALLEL_BOUNDARY || boundaryType == MOVING_WALL);

			//PARALLEL boundaries are not treated when they occur (e.g. at shared edges)
			for(int i = 0; i < Q && boundaryType != PARALLEL_BOUNDARY; ++i){
				int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
				int n_xyzoffset = p_computeNeighborCellOffset(k, j, fixedValue, c, xlength, wallIdx);
				int n_cell_index = Q*n_xyzoffset;

				//check (1) valid index: in case the direction of vector 'c' points to a non-existing cell
				//check (2) that the neighboring cell is a FLUID cell or PARALLEL_BOUNDARY (treated as inner)
				if(n_cell_index >= 0 && n_cell_index < maxValidIndex &&
						(flagField[n_xyzoffset] == FLUID || flagField[n_xyzoffset] == PARALLEL_BOUNDARY)
				){
					p_setBounceBack(collideField, procData.wallVelocity, boundaryType, i, current_cell_index, n_cell_index, c);
				}
			}
		}
	}
}

void treatBoundary(double * collideField, int const * const flagField, const t_procData procData){
	for(int wall=LEFT; wall<=BACK; ++wall){ //see LBDefinitions for order of walls
		if(procData.neighbours[wall] == MPI_PROC_NULL){
			// printWallEnum(wall);
			p_treatSingleWall(collideField, flagField, procData, wall);
		}
	}
}
