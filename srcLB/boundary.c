#include <stdio.h>
#include "boundary.h"
#include "helper.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

static inline int p_checkValidFluidIndex(const int totalSize, const int nextCellIndex,
	const int cellType) {
	return (nextCellIndex >= 0 && nextCellIndex < totalSize && cellType == FLUID);
}

void p_noSlip(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize){
	// Begin
	int i;
	int nextPoint[3];
	int currentCellIndex, nextFlagIndex, nextCellIndex;
	p_computeIndexQ(point, xlength, &currentCellIndex);
	for(i = 0; i < Q; i++){
		nextPoint[0] = point[0] + LATTICEVELOCITIES[i][0];
		nextPoint[1] = point[1] + LATTICEVELOCITIES[i][1];
		nextPoint[2] = point[2] + LATTICEVELOCITIES[i][2];
		p_computeIndex(nextPoint,  xlength, &nextFlagIndex);
		nextCellIndex = Q*nextFlagIndex;

	    if(p_checkValidFluidIndex(*totalSize, nextCellIndex, flagField[nextFlagIndex].type)) {
		    collideField[currentCellIndex + i] = collideField[nextCellIndex + (Q-i-1)];
		}
	}
}

void p_movingWall(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize) {
	// Begin
	int i;
	int nextPoint[3];
	int currentCellIndex, nextFlagIndex, nextCellIndex;
	double density;
	const double * const wallVelocity = boundPara->velocity;
	p_computeIndexQ(point, xlength, &currentCellIndex);

	for(i = 0; i < Q; i++){
		int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
		nextPoint[0] = point[0] + LATTICEVELOCITIES[i][0];
		nextPoint[1] = point[1] + LATTICEVELOCITIES[i][1];
		nextPoint[2] = point[2] + LATTICEVELOCITIES[i][2];
		p_computeIndex(nextPoint,  xlength, &nextFlagIndex);
		nextCellIndex = Q*nextFlagIndex;

	    if(p_checkValidFluidIndex(*totalSize, nextCellIndex, flagField[nextFlagIndex].type)) {
            double dot_uwall_c = wallVelocity[0]*c[0]+wallVelocity[1]*c[1]+wallVelocity[2]*c[2];
            double weight = LATTICEWEIGHTS[i];

            computeDensity(&collideField[nextCellIndex], &density);

            collideField[currentCellIndex + i] = collideField[nextCellIndex + (Q-i-1)] +
                    2 * weight * density * dot_uwall_c / (C_S*C_S);
        }
	}
}

// It is difficult since the indices are ordered lexicographically
// thereby destroying any kind of generalization

// Assumption: Only walls are allowed to be freeSlip. Obstacles are not allowed.

// Explanation: First we find the cell which is in the direction of the normal
// of the boundary cell. This depends on the wall to which the cell belongs.
// For example, for the YZ_BOTTOM or x = 0 wall - the normal is always (1,0,0)
// that is pointing in the +ve x direction.
// The variable nonZeroIndex stores the index for which the component of the normal
// is not zero or in other words is ±1. This is the component which will help us
// decide the mirrored index.
// This gives us the nextCell which shares a face with the boundaryCell in the
// direction of the normal.
// Once we have the index of the nextCell we check the flagField if it is a
// fluid cell.
// If it is a fluid cell then we simply loop over the distributions of the boundary
// cell using the index "i". For each index "i" we loop over the distributions again
// to find the index "j" which is the mirror image of "i" in the direction of the
// normal. This is achieved by checking if the "nonZeroIndex" component of the
// "i"th and the "j"th lattice velocity are opposites of each other "and" non zero
// "and" the other 2 components are equal. Once this index "j" is found we break
// the "j" loop and continue with the next "i". To reduce going over unneccesary "i"
// we check if the lattice velocity at "i" has a component in the direction of the
// normal. This is done using the projection and checking if it is +ve.

// void p_freeSlip(double* collideField, int const * const flagField,
// 	int const * const point, int const * const xlength,
// 	const t_boundPara * const boundPara, int const * const normal) {
// 	// Begin
// 	int i,j;
// 	int nextPoint[3] = {point[0]+normal[0],
// 						point[1]+normal[1],
// 						point[2]+normal[2]};
// 	int nonZeroIndex;
// 	if (abs(normal[0]) == 1) {
// 		nonZeroIndex = 0;
// 	} else if (abs(normal[1]) == 1) {
// 		nonZeroIndex = 1;
// 	} else if (abs(normal[2]) == 1) {
// 		nonZeroIndex = 2;
// 	}
//
// 	int currentFlagIndex, currentCellIndex, nextFlagIndex, nextCellIndex;
// 	p_computeIndex(point, xlength, &currentFlagIndex);
// 	p_computeIndex(nextPoint, xlength, &nextFlagIndex);
// 	currentCellIndex	= Q*currentFlagIndex;
// 	nextCellIndex 		= Q*nextFlagIndex;
// 	if (flagField[nextFlagIndex] == FLUID) {
// 		for (i = 0; i < Q; i++) {
// 			int c1[3] = {LATTICEVELOCITIES[i][0],
// 						 LATTICEVELOCITIES[i][1],
// 						 LATTICEVELOCITIES[i][2]};
// 			int projection = LATTICEVELOCITIES[i][0]*normal[0]+
// 							 LATTICEVELOCITIES[i][1]*normal[1]+
// 						 	 LATTICEVELOCITIES[i][2]*normal[2];
// 			if (projection >= 1) {
// 				for (j = 0; j < Q; j++) {
// 					int c2[3] = {LATTICEVELOCITIES[j][0],
// 								 LATTICEVELOCITIES[j][1],
// 								 LATTICEVELOCITIES[j][2]};
// 					if (c1[nonZeroIndex] == -c2[nonZeroIndex] &&
// 						c1[(nonZeroIndex+1)%3] == c2[(nonZeroIndex+1)%3] &&
// 						c1[(nonZeroIndex+2)%3] == c2[(nonZeroIndex+2)%3] &&
// 						c1[nonZeroIndex] == normal[nonZeroIndex]) {
// 						collideField[currentCellIndex+i] = collideField[nextCellIndex+j];
// 						break;
// 					}
// 				}
// 			}
// 		}
// 	}
// }

// MATLAB code used to generate the mirrored indices
// function generateMirrorIndices(idx,dir)
// % idx = 1 => x normal
// % idx = 2 => y normal
// % idx = 3 => z normal
// % dir = +1 => normal is in positive direction
// % dir = -1 => normal is in negative direction
// lv = [0,-1,-1;-1,0,-1;0,0,-1;1,0,-1;0,1,-1;
//    -1,-1,0;0,-1,0; 1,-1,0;-1,0,0;0,0,0;
//    1,0,0;  -1,1,0; 0,1,0; 1,1,0; 0,-1,1;
//    -1,0,1; 0,0,1;  1,0,1; 0,1,1];
// idx1 = idx;
// idx2 = mod(idx1,3)+1;
// idx3 = mod(idx1+1,3)+1;
// if idx2 == 0
//     idx2 = 1;
// end
// if idx3 == 0
//     idx3 = 1;
// end
// for i = 1:19
//     for j = 1:19
//         if lv(i,idx1) == -lv(j,idx1) && lv(i,idx1) == dir && lv(i,idx2) == lv(j,idx2) && lv(i,idx3) == lv(j,idx3)
//             disp([i-1,j-1])
//             break;
//         end
//     end
// end
// end

// Function to assign the indices
void p_assignIndices(const short int * const flag, int * index, int * mirrorIndex) {
	assert(*flag >= XY_LEFT && *flag <= XZ_BACK);

	switch (*flag) {
		case XY_LEFT:
			// z = 0
			index[0] = 14; index[1] = 15; index[2] = 16; index[3] = 17; index[4] = 18;
			mirrorIndex[0] = 0; mirrorIndex[1] = 1; mirrorIndex[2] = 2; mirrorIndex[3] = 3; mirrorIndex[4] = 4;
			break;
		case XY_RIGHT:
			// z = xlength[2]+1
			mirrorIndex[0] = 14; mirrorIndex[1] = 15; mirrorIndex[2] = 16; mirrorIndex[3] = 17; mirrorIndex[4] = 18;
			index[0] = 0; index[1] = 1; index[2] = 2; index[3] = 3; index[4] = 4;
			break;
		case YZ_BOTTOM:
			// x = 0
			index[0] = 3; index[1] = 7; index[2] = 10; index[3] = 13; index[4] = 17;
			mirrorIndex[0] = 1; mirrorIndex[1] = 5; mirrorIndex[2] = 8; mirrorIndex[3] =  11; mirrorIndex[4] = 15;
			break;
		case YZ_TOP:
			// x = xlength[0]+1
			mirrorIndex[0] = 3; mirrorIndex[1] = 7; mirrorIndex[2] = 10; mirrorIndex[3] = 13; mirrorIndex[4] = 17;
			index[0] = 1; index[1] = 5; index[2] = 8; index[3] =  11; index[4] = 15;
			break;
		case XZ_FRONT:
			// y = xlength[1]+1
			mirrorIndex[0] = 4; mirrorIndex[1] = 11; mirrorIndex[2] = 12; mirrorIndex[3] = 13; mirrorIndex[4] = 18;
			index[0] = 0; index[1] = 5; index[2] = 6; index[3] = 7; index[4] = 14;
			break;
		case XZ_BACK:
			// y = 0
			index[0] = 4; index[1] = 11; index[2] = 12; index[3] = 13; index[4] = 18;
			mirrorIndex[0] = 0; mirrorIndex[1] = 5; mirrorIndex[2] = 6; mirrorIndex[3] = 7; mirrorIndex[4] = 14;
			break;
		default:
			ERROR("** This should not happen!! **");
			break;
	}
}

void p_freeSlip(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize) {
	// Begin
	int i;
	int normal[3] = {0,0,0};
	int currentFlagIndex, currentCellIndex, nextFlagIndex, nextCellIndex;

	p_computeIndex(point, xlength, &currentFlagIndex);

	currentCellIndex	= Q*currentFlagIndex;

	// Assign normals based on the wall we are at
	/*TODO: (DL) readability bad... but saves the switch case. Can decide... once it works
	* I make a "speed difference test": */
//	static int normalIdx[6] = {2, 2, 0, 0 , 1, 1}; //static to allocate only once
//	static int normalVal[6] = {1,-1, 1,-1, -1, 1};
//	normal[ normalIdx[ flagField[currentFlagIndex].position ] ] =
//			normalVal[ flagField[currentFlagIndex].position ];

	switch (flagField[currentFlagIndex].position) {
		case XY_LEFT:
			normal[2] = 1;
			break;
		case XY_RIGHT:
			normal[2] = -1;
			break;
		case YZ_BOTTOM:
			normal[0] = 1;
			break;
		case YZ_TOP:
			normal[0] = -1;
			break;
		case XZ_BACK:
			normal[1] = 1;
			break;
		case XZ_FRONT:
			normal[1] = -1;
			break;
		default:
			ERROR("** This should not happen!! **");
			break;
	}

	int nextPoint[3] = {point[0]+normal[0],
						point[1]+normal[1],
						point[2]+normal[2]};

	p_computeIndex(nextPoint, xlength, &nextFlagIndex);
	nextCellIndex 		= Q*nextFlagIndex;

	int index[5] = {0,0,0,0,0}, mirrorIndex[5] = {0,0,0,0,0};

	p_assignIndices(&flagField[currentFlagIndex].position,index,mirrorIndex);

	//it's by construction always inside the domain, so no check for for loop needed
	assert(nextCellIndex >= 0 && nextCellIndex < *totalSize);

	//every cell is mirrored (not only FLUID)
	for (i = 0; i < 5; i++) {
		collideField[currentCellIndex+index[i]] = collideField[nextCellIndex+mirrorIndex[i]];
	}
}

void p_outflow(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize) {
	// Begin
	int i;
	int nextPoint[3];
	int nextFlagIndex, nextCellIndex, currentFlagIndex, currentCellIndex;
	double density, feq[Q], nextPointVel[3];

	p_computeIndex(point, xlength, &currentFlagIndex);
	currentCellIndex = Q*currentFlagIndex;

	for (i = 0; i < Q; i++) {
		nextPoint[0] = point[0] + LATTICEVELOCITIES[i][0];
		nextPoint[1] = point[1] + LATTICEVELOCITIES[i][1];
		nextPoint[2] = point[2] + LATTICEVELOCITIES[i][2];
		p_computeIndex(nextPoint, xlength, &nextFlagIndex);
		nextCellIndex = Q*nextFlagIndex;

		if(p_checkValidFluidIndex(*totalSize, nextCellIndex, flagField[nextFlagIndex].type)) {

			computeDensity(&collideField[nextCellIndex], &density);
			computeVelocity(&collideField[nextCellIndex], &density, nextPointVel);
			computeFeq(&(boundPara->rhoRef), nextPointVel, feq);
			collideField[currentCellIndex+i] = feq[Q-1-i] + feq[i] -
											collideField[nextCellIndex+Q-1-i];
		}
	}
}

void p_inflow(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize) {
	// Begin
	int i;
	int currentFlagIndex, currentCellIndex;
	double feq[Q];
	p_computeIndex(point, xlength, &currentFlagIndex);
	currentCellIndex = Q*currentFlagIndex;

	computeFeq(&boundPara->rhoRef, boundPara->velocity, feq);

	/* TODO: (DL) see WS sentence after eq. 2.2:
	 * "You may also try to use the density from the previous time step as pref"
	 *
	 * We could do this by either: using a static variable for pref or change the
	 * variable in boundPara (but then we need a non-const boundPara).
	 */

	for (i = 0; i < Q; i++) {
		collideField[currentCellIndex+i] = feq[i];
	}
}

void p_pressureIn(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize) {
	// Begin
	int i;
	int nextPoint[3];
	int nextFlagIndex, nextCellIndex, currentCellIndex;
	double feq[Q];

	p_computeIndexQ(point, xlength, &currentCellIndex);
	computeFeq(&boundPara->rhoIn, boundPara->velocity, feq);

	for (i = 0; i < Q; i++) {
		nextPoint[0] = point[0] + LATTICEVELOCITIES[i][0];
		nextPoint[1] = point[1] + LATTICEVELOCITIES[i][1];
		nextPoint[2] = point[2] + LATTICEVELOCITIES[i][2];
		p_computeIndex(nextPoint, xlength, &nextFlagIndex);
		nextCellIndex = Q*nextFlagIndex;

        if(p_checkValidFluidIndex(*totalSize, nextCellIndex, flagField[nextFlagIndex].type)) {
            collideField[currentCellIndex+i] = feq[Q-1-i] + feq[i] -
                                            collideField[nextCellIndex+Q-1-i];
        }
	}
}

t_boundaryFcnPtr p_selectFunction(const int wallType) {

	t_boundaryFcnPtr tmpPtr;
	switch (wallType) {
		case NO_SLIP:
			// call no slip wall
			tmpPtr = &p_noSlip;
			break;
		case MOVING:
			// call moving wall
			tmpPtr = &p_movingWall;
			break;
		case FREE_SLIP:
			// call free slip wall
			tmpPtr = &p_freeSlip;
			break;
		case OUTFLOW:
			// call outflow wall
			tmpPtr = &p_outflow;
			break;
		case INFLOW:
			// call inflow wall
			tmpPtr = &p_inflow;
			break;
		case PRESSURE_IN:
			// call pressure in wall
			tmpPtr = &p_pressureIn;
			break;
		case OBSTACLE:
			tmpPtr = &p_noSlip;
			break;
		default:
			ERROR("**** FLUID cell encountered! ****");
			tmpPtr = NULL;
			break;
	}
	assert(tmpPtr != NULL);
	return tmpPtr;
}


void treatBoundary(double *collideField, const t_flagField * const flagField,
	const t_boundPara * const boundPara, const int * const xlength){

	// iteration variables
	int x, y, z;
	int flagIndex, wallType, wallPos;
	int points[3];
	int const xlen2[3] = {xlength[0]+2,xlength[1]+2,xlength[2]+2};
    int const totalSize  = Q*xlen2[2]*xlen2[1]*xlen2[0];

	t_boundaryFcnPtr fcnPtr = NULL;

	// TODO: (VS) Compute obstacle min and max to exclude some cells
	// As of now iterate over everything normally

	// TODO: (VS) directly call functions rather than having indirect call

	for (z = 0; z < xlen2[2]; z++) {
		for (y = 0; y < xlen2[1]; y++) {
			for (x = 0; x < xlen2[0]; x++) {
				p_computeIndexXYZ(x,y,z,xlength,&flagIndex);
				wallType	= flagField[flagIndex].type;
				wallPos		= flagField[flagIndex].position;

				assert(wallType != INVALID);
				if (wallType != FLUID) {
					points[0] = x;
					points[1] = y;
					points[2] = z;
					fcnPtr = p_selectFunction(wallType);
					(*fcnPtr)(collideField, flagField, points, xlength,
						&boundPara[wallPos], &totalSize);
				}
			}
		}
	}
}
