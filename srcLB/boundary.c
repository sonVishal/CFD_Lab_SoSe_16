#include <stdio.h>
#include "boundary.h"
#include "helper.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

// Function that checks the index is valid and if it is valid then
// checks whether the cell is FLUID or not
static inline int p_checkValidFluidIndex(const int totalSize, const int nextCellIndex, const t_flagField * const flagField) {
	if (nextCellIndex >= 0 && nextCellIndex < totalSize) {
		if (flagField[nextCellIndex].type == FLUID) {
			return 1;
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}

// Handle no slip boundary condition
void p_noSlip(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize){
	// Begin
	int i;
	int nextPoint[3];
	int currentCellIndex, nextFlagIndex, nextCellIndex;
	// Compute the index of the current cell
	p_computeIndexQ(point, xlength, &currentCellIndex);

	// Loop over the distributions
	for(i = 0; i < Q; i++){
		// Compute the next cell coordinates
		nextPoint[0] = point[0] + LATTICEVELOCITIES[i][0];
		nextPoint[1] = point[1] + LATTICEVELOCITIES[i][1];
		nextPoint[2] = point[2] + LATTICEVELOCITIES[i][2];
		p_computeIndex(nextPoint,  xlength, &nextFlagIndex);
		nextCellIndex = Q*nextFlagIndex;
		// Check if fluid and assign the inverse distribution from next cell
	    if(p_checkValidFluidIndex(*totalSize, nextCellIndex, flagField)) {
		    collideField[currentCellIndex + i] = collideField[nextCellIndex + (Q-i-1)];
		}
	}
}

// Handle moving wall boundary condition
void p_movingWall(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize) {
	// Begin
	int i;
	int nextPoint[3];
	int currentCellIndex, nextFlagIndex, nextCellIndex;
	double density;
	// Store the wall velocity
	const double * const wallVelocity = boundPara->velocity;
	// Compute the current cell index
	p_computeIndexQ(point, xlength, &currentCellIndex);

	// Loop over all the distribution
	for(i = 0; i < Q; i++){
		int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
		// Compute next cell coordinates
		nextPoint[0] = point[0] + LATTICEVELOCITIES[i][0];
		nextPoint[1] = point[1] + LATTICEVELOCITIES[i][1];
		nextPoint[2] = point[2] + LATTICEVELOCITIES[i][2];
		// Compute next cell index
		p_computeIndex(nextPoint,  xlength, &nextFlagIndex);
		nextCellIndex = Q*nextFlagIndex;

		// Check if fluid and assign the wall velocity
	    if(p_checkValidFluidIndex(*totalSize, nextCellIndex, flagField)) {
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

// Function to assign the indices and mirror indices
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
	}
}

// Handle free slip boundary condition
void p_freeSlip(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize) {
	// Begin
	int i;
	int normal[3] = {0,0,0};
	int currentFlagIndex, currentCellIndex, nextFlagIndex, nextCellIndex;
	// Compute current cell index
	p_computeIndex(point, xlength, &currentFlagIndex);
	currentCellIndex	= Q*currentFlagIndex;

	// Assign normals based on the wall we are at
	// switch (flagField[currentFlagIndex].position) {
	// 	case XY_LEFT:
	// 		normal[2] = 1;
	// 		break;
	// 	case XY_RIGHT:
	// 		normal[2] = -1;
	// 		break;
	// 	case YZ_BOTTOM:
	// 		normal[0] = 1;
	// 		break;
	// 	case YZ_TOP:
	// 		normal[0] = -1;
	// 		break;
	// 	case XZ_BACK:
	// 		normal[1] = 1;
	// 		break;
	// 	case XZ_FRONT:
	// 		normal[1] = -1;
	// 		break;
	// 	default:
	// 		ERROR("** This should not happen!! **");
	// 		break;
	// }

	// Assert that the position is a wall
	assert(flagField[currentFlagIndex].position >= XY_LEFT &&
			flagField[currentFlagIndex].position <= XZ_BACK);

	// The index for the direction of the normal
	// static to allocate only once
	static int normalIdx[6] = {2, 2, 0, 0 , 1, 1};
	// The direction of the normal
	static int normalVal[6] = {1,-1, 1,-1, -1, 1};
	// Assign the respective normal
	normal[normalIdx[flagField[currentFlagIndex].position]] =
			normalVal[flagField[currentFlagIndex].position];

	// Compute the next point
	int nextPoint[3] = {point[0]+normal[0],
						point[1]+normal[1],
						point[2]+normal[2]};

	// Compute next cell index
	p_computeIndex(nextPoint, xlength, &nextFlagIndex);
	nextCellIndex 		= Q*nextFlagIndex;

	// Allocate memory
	int index[5] = {0,0,0,0,0}, mirrorIndex[5] = {0,0,0,0,0};

	// Get the index and mirror index for the distribution
	p_assignIndices(&flagField[currentFlagIndex].position,index,mirrorIndex);

	// it's by construction always inside the domain, so no check for for-loop needed
	assert(nextCellIndex >= 0 && nextCellIndex < Q*(*totalSize));

	// every cell is mirrored (not only FLUID)
	// Loop over only the 5 distributions
	// We need to do this because there are only 5 distributions per face
	for (i = 0; i < 5; i++) {
		collideField[currentCellIndex+index[i]] = collideField[nextCellIndex+mirrorIndex[i]];
	}
}

// Handle outflow boundary condition
void p_outflow(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize) {
	// Begin
	int i;
	int nextPoint[3];
	int nextFlagIndex, nextCellIndex, currentFlagIndex, currentCellIndex;
	double density, feq[Q], nextCellVel[3];

	// Compute current cell index
	p_computeIndex(point, xlength, &currentFlagIndex);
	currentCellIndex = Q*currentFlagIndex;

	// Loop over all the distribution
	for (i = 0; i < Q; i++) {
		// Compute next cell index
		nextPoint[0] = point[0] + LATTICEVELOCITIES[i][0];
		nextPoint[1] = point[1] + LATTICEVELOCITIES[i][1];
		nextPoint[2] = point[2] + LATTICEVELOCITIES[i][2];
		p_computeIndex(nextPoint, xlength, &nextFlagIndex);
		nextCellIndex = Q*nextFlagIndex;

		// If valid fluid cell then apply the boundary condition
		if(p_checkValidFluidIndex(*totalSize, nextCellIndex, flagField)) {

			computeDensity(&collideField[nextCellIndex], &density);
			computeVelocity(&collideField[nextCellIndex], &density, nextCellVel);
			computeFeq(&(boundPara->rhoRef), nextCellVel, feq);
			collideField[currentCellIndex+i] = feq[Q-1-i] + feq[i] -
											collideField[nextCellIndex+Q-1-i];
		}
	}
}

// Handle the inflow boundary condition
void p_inflow(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize) {
	// Begin
	// Counter so that inflow is called only once
	static int counter = 0;
	if (counter == 0) {
		int i;
		int currentFlagIndex, currentCellIndex;
		double feq[Q];
		// Compute the current cell index
		p_computeIndex(point, xlength, &currentFlagIndex);
		currentCellIndex = Q*currentFlagIndex;

		// Compute equilibrium distribution
		computeFeq(&boundPara->rhoRef, boundPara->velocity, feq);

		// Assign the distributions
		for (i = 0; i < Q; i++) {
			collideField[currentCellIndex+i] = feq[i];
		}

		// Increment counter to avoid calling this function again
		counter++;
	} else {
		return;
	}
}

// Handle the pressure in boundary condition
void p_pressureIn(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize) {
	// Begin
	int i;
	int nextPoint[3];
	int nextFlagIndex, nextCellIndex, currentFlagIndex, currentCellIndex;
	double density, feq[Q], nextCellVel[3];
	// Compute the current cell index
	p_computeIndex(point, xlength, &currentFlagIndex);
	currentCellIndex = Q*currentFlagIndex;

	// Loop over all the distributions
	for (i = 0; i < Q; i++) {
		// Compute next cell index
		nextPoint[0] = point[0] + LATTICEVELOCITIES[i][0];
		nextPoint[1] = point[1] + LATTICEVELOCITIES[i][1];
		nextPoint[2] = point[2] + LATTICEVELOCITIES[i][2];
		p_computeIndex(nextPoint, xlength, &nextFlagIndex);
		nextCellIndex = Q*nextFlagIndex;

		// If fluid then apply the boundary condition
		if(p_checkValidFluidIndex(*totalSize, nextCellIndex, flagField)) {

			computeDensity(&collideField[nextCellIndex], &density);
			computeVelocity(&collideField[nextCellIndex], &density, nextCellVel);
			computeFeq(&(boundPara->rhoIn), nextCellVel, feq);
			collideField[currentCellIndex+i] = feq[Q-1-i] + feq[i] -
											collideField[nextCellIndex+Q-1-i];
		}
	}
}

// Choose the function based on the boundary type
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

// Treat the boundary
void treatBoundary(double *collideField, const t_flagField * const flagField,
	const t_boundPara * const boundPara, const int * const xlength){
	// Begin
	// iteration variables
	int x, y, z;
	int flagIndex, wallType, wallPos;
	int points[3];
	int const xlen2[3] = {xlength[0]+2,xlength[1]+2,xlength[2]+2};
    int const totalSize  = xlen2[2]*xlen2[1]*xlen2[0];

	t_boundaryFcnPtr fcnPtr = NULL;

	// Iterate over everything
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
