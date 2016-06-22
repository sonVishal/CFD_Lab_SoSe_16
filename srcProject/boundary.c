#include "boundary.h"

//TODO: (DL) in case we do not need the boundaries NO_SLIP and MOVING_WALL then delete the code for clarity.
//However, leave it as long as possible for testing and compare with previous scenarios.

/* Helper function that carries out the moving wall or bounce back condition, depending on flag 'type'.
 * The value in flagField corresponding to the current cell (current_cell_index) has to indicate a boundary
 * cell (flag 1 or 2), the value in flagField corresponding to neighboring cell (n_cell_index) has to
 * be a cell in the domain (= FLUID cell of flag 0).
 */
//void p_setBounceBack(double *collideField, const double * const wallVelocity,
		//const int type, const int i, const int current_cell_index, const int n_cell_index, const int *c) {

	//if(type == 1){ 	 //bounce back
		//[> equation 17 <]
		//collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)];
	//}else{ // type == 2 (MOVING_WALL)

		//double density = 0;
		//computeDensity(&collideField[n_cell_index], &density);

		//double dot_uwall_c = wallVelocity[0]*c[0]+wallVelocity[1]*c[1]+wallVelocity[2]*c[2];
		//double weight = LATTICEWEIGHTS[i];

		//[> equation 19 <]
		//collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)] +
				//2 * weight * density * dot_uwall_c / (C_S*C_S);
	//}
//}

/* Helper function similar to p_computeCellOffset, except that it computes the offset value
 * of a neighboring cell in the direction of 'c_vec'.
 * There are no checks that guarantee that the computed offset value is valid.
 * For example, when the current cell is on the boundary and 'c_vec' points away from the inner value.
 */
//int p_computeNeighborCellOffset(int outer, int inner, int fixedValue,
		//const int * const c_vec, int const * const xlength, const int wallIdx){

	//switch (wallIdx/2) {
	//case 0: // LEFT, RIGHT -> Y fixed
		//outer 		+= c_vec[2];	//= Z
		//inner 		+= c_vec[0];	//= X
		//fixedValue 	+= c_vec[1];	//= Y
		//return (xlength[0]+2) * (outer * (xlength[1]+2) + fixedValue) + inner;

	//case 1: // TOP, BOTTOM -> Z fixed
		//outer 		+= c_vec[1];	//= Y
		//inner 		+= c_vec[0];	//= X
		//fixedValue 	+= c_vec[2];	//= Z
		//return (xlength[0]+2) * (fixedValue * (xlength[1]+2) + outer) + inner;

	//case 2: // FRONT, BACK -> X fixed
		//outer 		+= c_vec[2];	//= Z
		//inner 		+= c_vec[1];	//= Y
		//fixedValue 	+= c_vec[0];	//= X
		//return (xlength[0]+2) * (outer * (xlength[1]+2) + inner) + fixedValue;

	//default:
		//ERROR("Invalid wall index occured. This should not happen !!!");
		//return -1;
	//}
//}


/* Helper function that treats a single boundary wall. The wall itself is described with wallIdx
 * (from enum WALLS) and the 'fixedValue' (generally 0 or xlength+1).
 */
//void p_treatSingleWall(double *collideField, const int * const flagField, const t_procData * const procData, const int wallIdx){

	//int endOuter = -1, endInner = -1, fixedValue = -1;
	//p_setIterationParameters(&endOuter, &endInner, &fixedValue, procData, wallIdx);
	//assert(endOuter != -1 && endInner != -1 && fixedValue != -1);

	//int xlength[3] = {procData->xLength[0], procData->xLength[1], procData->xLength[2]};

	////variable needed at check whether it is an in-domain (FLUID) cell
	//int maxValidIndex = Q*(xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2);

	////k - corresponds to the 'outer' value when computing the offset
	//for(int k = 0; k <= endOuter; ++k){
		////j - corresponds to the 'inner' value
		//for(int j = 0; j <= endInner; ++j){

			////flagField[xyz_offset] should only be a boundary cell (checked in p_setBounceBack)
			//int xyz_offset = p_computeCellOffset(k, j, fixedValue, xlength, wallIdx);
			//int current_cell_index = Q*xyz_offset;
			//int boundaryType = flagField[xyz_offset];

			//assert(boundaryType == NO_SLIP || boundaryType == MOVING_WALL);

			////PARALLEL boundaries are not treated when they occur (e.g. at shared edges)
			//for(int i = 0; i < Q; ++i){
				//int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
				//int n_xyzoffset = p_computeNeighborCellOffset(k, j, fixedValue, c, xlength, wallIdx);
				//int n_cell_index = Q*n_xyzoffset;

				////check (1) valid index: in case the direction of vector 'c' points to a non-existing cell
				////check (2) that the neighboring cell is a FLUID cell or PARALLEL_BOUNDARY (treated as inner)
				//if(n_cell_index >= 0 && n_cell_index < maxValidIndex &&
						//(flagField[n_xyzoffset] == FLUID || flagField[n_xyzoffset] == PARALLEL_BOUNDARY)
				//){
					//p_setBounceBack(collideField, procData->wallVelocity, boundaryType, i, current_cell_index, n_cell_index, c);
				//}
			//}
		//}
	//}
//}
#define VARIABLE -1 //symbol that is used in context of indices. It indicates which index is variable.

int extractInjectEdge(double buffer[], double * const collideField, t_iterParaEdge const * const iterPara,
		t_procData const * const procData, const int index, const int injectFlag){

	assert(injectFlag == 0 || injectFlag == 1);

	int endIndex, cellOffset;

	//Note: the corners of shared edges do not have to be treated
	//(in D3Q19 there are no distributions into the corners)
	if(iterPara->x == VARIABLE){
		endIndex = procData->xLength[0];
	}else if(iterPara->y == VARIABLE){
		endIndex = procData->xLength[1];
	}else{
		assert(iterPara->z == VARIABLE);
		endIndex = procData->xLength[2];
	}

	int currentIndexBuff = 0, currentIndexField = -1;

	//z * (xlen*ylen) + y * (xlen) + x
	for(int varIdx = 1; varIdx <= endIndex; ++varIdx){
		if(iterPara->x == VARIABLE){
			cellOffset = procData->xLength[0] * (procData->xLength[1] * iterPara->z + iterPara->y) + varIdx;
		}else if(iterPara->y == VARIABLE){
			cellOffset = procData->xLength[0] * (procData->xLength[1] * iterPara->z + varIdx) + iterPara->x;
		}else{
			assert(iterPara->z == VARIABLE);
			cellOffset = procData->xLength[0] * (procData->xLength[1] * varIdx + iterPara->y) + iterPara->x;
		}

		currentIndexField = Q*cellOffset;
		if(injectFlag)
			collideField[currentIndexField+index] = buffer[currentIndexBuff++];
		else
			buffer[currentIndexBuff++] = collideField[currentIndexField+index];
	}
	return currentIndexBuff-1; //Will be used as bufferSize
}


void p_setBoundaryIterParameters(t_iterPara *const iterPara, t_procData const*const procData, const int direction){

	switch(direction/2){ //integer division to get the type of face (see enum in LBDefinitions.h)
	iterPara->startOuter = 1;
	iterPara->startInner = 1;

	//---------------------------------------------
	//outer = Z, inner = X, Y fixed
	case 0:
		iterPara->endOuter   = procData->xLength[2];
		iterPara->endInner   = procData->xLength[0];
		iterPara->fixedValue = (direction == LEFT) ? 1 : procData->xLength[1];
		break;

	//---------------------------------------------
	//outer = Y, inner = X, Z fixed
	case 1:
		iterPara->endOuter   = procData->xLength[1];
		iterPara->endInner   = procData->xLength[0];
		iterPara->fixedValue = (direction == BOTTOM) ? 1 : procData->xLength[2];
		break;

	//---------------------------------------------
	//outer = Z, inner = Y, X fixed
	case 2:
		iterPara->endOuter   = procData->xLength[2];
		iterPara->endInner   = procData->xLength[1];
		iterPara->fixedValue = (direction == BACK) ? 1 : procData->xLength[0];
		break;

	default:
		ERROR("Invalid direction occurred. This should never happen!");
	}
}

void p_setEdgeIterParameters(t_iterParaEdge *const iterPara, t_procData const*const procData, const int edge, const int injectFlag){

	assert(injectFlag == 0 || injectFlag == 1);

	//if inject: copy into boundary; else: copy from domain cells lying at edge
	int start     = injectFlag ? 0 : 1; //if inject: start = 0, is at (ghost) boundary layer
	int endOffset = injectFlag ? 1 : 0; //if inject: add endOffset +1 to be at (ghost) boundary layer

	switch(edge/4){ //integer division - 12 edges, 3 cases
	case 0: //edges 0-3 (horizontal, bottom)
		iterPara->z = 0;
		if(edge == 0 || edge == 2){
			iterPara->x = VARIABLE;
			iterPara->y = (edge == 0) ? start : procData->xLength[1]+endOffset;
		}else{ //edge == 1 || edge == 3
			iterPara->y = VARIABLE;
			iterPara->x = (edge == 1) ? procData->xLength[0]+endOffset : start;
		}
		break;

	case 1: //edges 4-7 (vertical)
		iterPara->z = VARIABLE;
		iterPara->y = (edge == 4 || edge == 5) ? start : procData->xLength[1]+endOffset;
		iterPara->x = (edge == 4 || edge == 7) ? start : procData->xLength[0]+endOffset;
		break;

	case 2: //edges 8-11 (horizontal, top)
		iterPara->z = procData->xLength[2]+endOffset;
		if(edge == 8 || edge == 10){
			iterPara->x = VARIABLE;
			iterPara->y = (edge == 8) ? start : procData->xLength[1]+endOffset;
		}else{ //edge == 9 || edge == 11
			iterPara->y = VARIABLE;
			iterPara->x = (edge == 9) ? procData->xLength[0]+endOffset : start;
		}
		break;

	default:
		assert(edge>=0 && edge<=11);
	}
	assert(((iterPara->x < 0) + (iterPara->y < 0) + (iterPara->z < 0)) == 1); //only one is VARIABLE
}

int p_assignSharedEdgeIndex(const int edge) {

	//Look for edges numbering in LBDefinitions.h
	// edge 0: (0,-1,-1) 	-> [0]
	// edge 1: (1,0,-1) 	-> [3]
	// edge 2: (0,1,-1) 	-> [4]
	// edge 3: (-1,0,-1) 	-> [1]

	// edge 4: (-1,-1,0) 	-> [5]
	// edge 5: (1,-1,0) 	-> [7]
	// edge 6: (1,1,0) 		-> [13]
	// edge 7: (-1,1,0) 	-> [11]

	// edge 8: (0,-1,1) 	-> [14]
	// edge 9: (1,0,1) 		-> [17]
	// edge 10: (0,1,1) 	-> [18]
	// edge 11: (-1,0,1) 	-> [15]
	static int indices[12] = {0,3,4,1,
							  5,7,13,11,
							  14,17,18,15};
	return indices[edge];
}


void treatPeriodicWall(double *const collideField, double *const sendBuffer, double *const readBuffer,
	const t_procData * const procData, const int procWall, const int opponentWall){

	t_iterPara  iterPara;
	int 		indexIn[5], indexOut[5], bufferSize, commRank;

	p_assignIndices(procWall,     indexOut);
	p_assignIndices(opponentWall, indexIn);

	p_setBoundaryIterParameters(&iterPara, procData, procWall);

	//TODO: (DL) possibly safe this somewhere...
	//always excluding shared edges (treat differently&separately)
	bufferSize = 5 * iterPara.endOuter * iterPara.endInner;

	//TODO: extract and inject are from parallel.c - Discuss if we make them "common functions"
	extract(sendBuffer, collideField, &iterPara, procData, procWall, indexOut);

	commRank = procData->periodicNeighbours[procWall];

	//tag is chosen to be 1 that it is different from parallel_boundary (0) and edges (2)
	MPI_Sendrecv(sendBuffer, bufferSize, MPI_DOUBLE, commRank, 1, readBuffer,
		bufferSize, MPI_DOUBLE, commRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	inject(readBuffer, collideField, &iterPara, procData, procWall, indexIn);
}

void treatPeriodicEdge(double *collideField, double *const sendBuffer, double *const readBuffer,
	const t_procData * const procData, const int procEdge, const int opponentEdge){

	t_iterParaEdge iterParaEdge;
	//Note: in edge case only one index has to be copied
	int indexOutEdge, indexInEdge, bufferSize, commRank;

	indexOutEdge = p_assignSharedEdgeIndex(procEdge);
	indexInEdge =  p_assignSharedEdgeIndex(opponentEdge);

	p_setEdgeIterParameters(&iterParaEdge, procData, procEdge, 0); //0 = extract

	//using buffers from parallel boundaries
	bufferSize = extractInjectEdge(sendBuffer, collideField, &iterParaEdge, procData, indexOutEdge, 0);

	commRank = procData->periodicEdgeNeighbours[procEdge]; //communication rank

	//tag is chosen to be 2 that it is different from parallel_boundary (0) and periodic_boundary (1)
	MPI_Sendrecv(sendBuffer, bufferSize, MPI_DOUBLE, commRank, 2, readBuffer,
			bufferSize, MPI_DOUBLE, commRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	p_setEdgeIterParameters(&iterParaEdge, procData, procEdge, 1); //1 = inject
	extractInjectEdge(readBuffer, collideField, &iterParaEdge, procData, indexInEdge, 1);
}

void treatPeriodicWallNoComm(){
	ERROR("implement");
}

void treatPeriodicEdgeNoComm(){
	ERROR("implement");
}


void treatComponentBoundary(t_component *c,int numComp, int const * const flagField, const t_procData * const procData, double **sendBuffer, double **readBuffer){
    for (int i = 0; i < numComp; ++i) {
        treatBoundary(flagField, c[i].collideField, procData, sendBuffer, readBuffer);
    }
}

//TODO: (TKS) Remove flagField if not in use when finished
//TODO: (DL) Update header with new functions if necessary.
void treatBoundary(int const * const flagField, double *collideField, const t_procData * const procData, double **sendBuffer, double **readBuffer){

	// Handle of "real" boundaries:
	// const int NO_NEIGHBOUR = -2; // equals the MPI_PROC_NULL = -2
	// for(int wall=LEFT; wall<=BACK; ++wall){ //see LBDefinitions for order of walls
	// 	if(procData->neighbours[wall] == NO_NEIGHBOUR){
	// 		// printWallEnum(wall);
	// 		p_treatSingleWall(c->collideField, flagField, procData, wall);
	// 	}
	// }

	for(int wall=LEFT; wall<=BACK; wall+=2){ //see LBDefinitions for order of walls

		int firstNeighbour = procData->periodicNeighbours[wall];
		int secondNeighbour = procData->periodicNeighbours[wall+1];

		if(firstNeighbour != MPI_PROC_NULL && secondNeighbour != MPI_PROC_NULL){
			//Case when there is only one proc, then no communication is required
			treatPeriodicWallNoComm(); //TODO: (DL) implement
		}
		else{
			if(firstNeighbour != MPI_PROC_NULL){
				treatPeriodicWall(collideField, sendBuffer[wall], readBuffer[wall],
					procData, wall, wall+1);
			}

			if(secondNeighbour != MPI_PROC_NULL){
				treatPeriodicWall(collideField, sendBuffer[wall], readBuffer[wall],
					procData, wall+1, wall);
			}
		}
	}

	static const int edge1[6] = {0,1,2,3,4,5};
	static const int edge2[6] = {10,11,8,9,6,7}; //opposite edge to edge1; see numbering in LBDefinitions.h

	for(int idx = 0; idx < 6; ++idx){

		int firstNeighbour = procData->periodicEdgeNeighbours[edge1[idx]];
		int secondNeighbour = procData->periodicEdgeNeighbours[edge2[idx]];

		if(firstNeighbour != MPI_PROC_NULL && secondNeighbour != MPI_PROC_NULL){
			//Case when there is only one proc, then no communication is required
			treatPeriodicEdgeNoComm(); //TODO: (DL) implement
		}else{
			if(firstNeighbour != MPI_PROC_NULL){
				treatPeriodicEdge(collideField, sendBuffer[idx], readBuffer[idx],
					procData, edge1[idx], edge2[idx]);
			}

			if(secondNeighbour != MPI_PROC_NULL){
				treatPeriodicEdge(collideField, sendBuffer[idx], readBuffer[idx],
					procData, edge2[idx], edge1[idx]);
			}
		}
	}
}
