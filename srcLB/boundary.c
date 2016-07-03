#include "boundary.h"

/* Enum that gives a distinction of the walls with respective to the fixed value of (x,y,z) */
enum WALLS{E_XFIXED, E_YFIXED, E_ZFIXED};

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
	}else if(type == 2){ 	//moving wall

		double density = 0;
		computeNumDensity(&collideField[n_cell_index], &density);

		double dot_uwall_c = wallVelocity[0]*c[0]+wallVelocity[1]*c[1]+wallVelocity[2]*c[2];
		double weight = LATTICEWEIGHTS[i];

		/* equation 19 */
		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)] +
				2 * weight * dot_uwall_c / (C_S*C_S);

	}else{ 		//Wrong input parameter
		ERROR("A FLUID cell appeared when setting boundaries. This should not happen!!\n");
	}
}

/* Helper function that computes the offset of the current cell. 'Inner' corresponds to the first value of (x,y,z)
 * that is not fixed; 'outer' to the second value. By this ordering a better cache efficiency is obtained.
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
			outer 		+= c_vec[2];	//= Z
			inner 		+= c_vec[1];	//= Y
			fixedValue 	+= c_vec[0];	//= X
			return (xlength+2) * (outer*(xlength+2) + inner) + fixedValue;
		case E_YFIXED:
			outer 		+= c_vec[2];	//= Z
			inner 		+= c_vec[0];	//= X
			fixedValue 	+= c_vec[1];	//= Y
			return (xlength+2) * (outer*(xlength+2) + fixedValue) + inner;
		case E_ZFIXED:
			outer 		+= c_vec[1];	//= Y
			inner 		+= c_vec[0];	//= X
			fixedValue 	+= c_vec[2];	//= Z
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

	//variable needed at check whether it is an in-domain (FLUID) cell
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

				//check (1) valid index: in case the direction of vector 'c' points to a non-existing cell
				//check (2) that the neighboring cell is a FLUID cell in the domain (and not another boundary cell)
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

void treatBoundary(t_component *c, int* flagField, const double * const wallVelocity, int xlength, int includingDensity){
	for (int i = 0; i < NUMCOMP; i++) {
		treatWallPeriodic(c[i].collideField, c[i].rho, LEFT, xlength, includingDensity);
		treatWallPeriodic(c[i].collideField, c[i].rho, RIGHT, xlength, includingDensity);
		treatWallPeriodic(c[i].collideField, c[i].rho, TOP, xlength, includingDensity);
		treatWallPeriodic(c[i].collideField, c[i].rho, BOTTOM, xlength, includingDensity);
		treatWallPeriodic(c[i].collideField, c[i].rho, FRONT, xlength, includingDensity);
		treatWallPeriodic(c[i].collideField, c[i].rho, BACK, xlength, includingDensity);
	}
}

// This will be useful for computing otherSideIdx
// nbhR = (myRank+1)%numRanks;
// nbhL = (myRank-1+numRanks)%numRanks;
void treatWallPeriodic(double * collideField, double *rho, int direction, int xlength, int includingDensity) {
	int startOuter=0, startInner=0, endOuter=0, endInner=0, fixedValueIdx=0, fixedValueOtherIdx=0;
	switch (direction) {
		case LEFT:
			// z
			startOuter = 1; endOuter = xlength;
			// x
			startInner = 1; endInner = xlength;
			// y
			fixedValueIdx = 1;
			fixedValueOtherIdx = xlength+1;
			break;
		case RIGHT:
			// z
			startOuter = 1; endOuter = xlength;
			// x
			startInner = 1; endInner = xlength;
			// y
			fixedValueIdx = xlength;
			fixedValueOtherIdx = 0;
			break;
		case TOP:
			// y
			startOuter = 0; endOuter = xlength+1;
			// x
			startInner = 1; endInner = xlength;
			// z
			fixedValueIdx = xlength;
			fixedValueOtherIdx = 0;
			break;
		case BOTTOM:
			// y
			startOuter = 0; endOuter = xlength+1;
			// x
			startInner = 1; endInner = xlength;
			// z
			fixedValueIdx = 1;
			fixedValueOtherIdx = xlength+1;
			break;
		case FRONT:
			// z
			startOuter = 0; endOuter = xlength+1;
			// y
			startInner = 0; endInner = xlength+1;
			// x
			fixedValueIdx = xlength;
			fixedValueOtherIdx = 0;
			break;
		case BACK:
			// z
			startOuter = 0; endOuter = xlength+1;
			// y
			startInner = 0; endInner = xlength+1;
			// x
			fixedValueIdx = 1;
			fixedValueOtherIdx = xlength+1;
			break;
		default:
			ERROR("NO!");
			break;
	}
	for (int k = startOuter; k <= endOuter; k++) {
		for (int j = startInner; j <= endInner; j++) {
			int idx = Q*computeCellOffset(k, j, fixedValueIdx, direction, xlength);
			int otherSideIdx = Q*computeCellOffset(k, j, fixedValueOtherIdx, direction, xlength);

			if(includingDensity){
				rho[otherSideIdx/Q] = rho[idx/Q];
			}

			for (int i = 0; i < Q; i++) {
				collideField[otherSideIdx+i] = collideField[idx+i];
			}
		}
	}
}

int computeCellOffset(const int outer, const int inner, const int fixed, const int dir, const int xlength) {
	switch (dir) {
		case LEFT:
			return p_computeCellOffsetXYZ(inner, fixed, outer, xlength);
		case RIGHT:
			return p_computeCellOffsetXYZ(inner, fixed, outer, xlength);
		case TOP:
			return p_computeCellOffsetXYZ(inner, outer,fixed, xlength);
		case BOTTOM:
			return p_computeCellOffsetXYZ(inner, outer,fixed, xlength);
		case FRONT:
			return p_computeCellOffsetXYZ(fixed, inner, outer,xlength);
		case BACK:
			return p_computeCellOffsetXYZ(fixed, inner, outer,xlength);
		default:
			return -1;
			ERROR("NO!");
			break;
	}
}
