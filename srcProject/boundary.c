	#include "boundary.h"

//TODO: (DL) in case we do not need the boundaries NO_SLIP and MOVING_WALL then delete the code for clarity.
//However, leave it as long as possible for testing and compare with previous scenarios.

//TODO: (DL) reminder: delete all the printfs for debugging when finalizing

/* Helper function that carries out the moving wall or bounce back condition, depending on flag 'type'.
 * The value in flagField corresponding to the current cell (current_cell_index) has to indicate a boundary
 * cell (flag 1 or 2), the value in flagField corresponding to neighboring cell (n_cell_index) has to
 * be a cell in the domain (= FLUID cell of flag 0).
 */
// void p_setBounceBack(double *collideField, const double * const wallVelocity,
// 		const int type, const int i, const int current_cell_index, const int n_cell_index, const int *c) {
//
// 	if(type == 1){ 	 //bounce back
// 		// [> equation 17 <]
// 		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)];
// 	}else{ // type == 2 (MOVING_WALL)
//
// 		double density = 0;
// 		computeDensity(&collideField[n_cell_index], &density);
//
// 		double dot_uwall_c = wallVelocity[0]*c[0]+wallVelocity[1]*c[1]+wallVelocity[2]*c[2];
// 		double weight = LATTICEWEIGHTS[i];
//
// 		// [> equation 19 <]
// 		collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)] +
// 				2 * weight * density * dot_uwall_c / (C_S*C_S);
// 	}
// }

/* Helper function similar to p_computeCellOffset, except that it computes the offset value
 * of a neighboring cell in the direction of 'c_vec'.
 * There are no checks that guarantee that the computed offset value is valid.
 * For example, when the current cell is on the boundary and 'c_vec' points away from the inner value.
 */
// int p_computeNeighborCellOffset(int outer, int inner, int fixedValue,
// 		const int * const c_vec, int const * const xlength, const int wallIdx){
//
// 	switch (wallIdx/2) {
// 	case 0: // LEFT, RIGHT -> Y fixed
// 		outer 		+= c_vec[2];	//= Z
// 		inner 		+= c_vec[0];	//= X
// 		fixedValue 	+= c_vec[1];	//= Y
// 		return (xlength[0]+2) * (outer * (xlength[1]+2) + fixedValue) + inner;
//
// 	case 1: // TOP, BOTTOM -> Z fixed
// 		outer 		+= c_vec[1];	//= Y
// 		inner 		+= c_vec[0];	//= X
// 		fixedValue 	+= c_vec[2];	//= Z
// 		return (xlength[0]+2) * (fixedValue * (xlength[1]+2) + outer) + inner;
//
// 	case 2: // FRONT, BACK -> X fixed
// 		outer 		+= c_vec[2];	//= Z
// 		inner 		+= c_vec[1];	//= Y
// 		fixedValue 	+= c_vec[0];	//= X
// 		return (xlength[0]+2) * (outer * (xlength[1]+2) + inner) + fixedValue;
//
// 	default:
// 		ERROR("Invalid wall index occured. This should not happen !!!");
// 		return -1;
// 	}
// }

/* Helper function that treats a single boundary wall. The wall itself is described with wallIdx
 * (from enum WALLS) and the 'fixedValue' (generally 0 or xlength+1).
 */
// void p_treatSingleWall(double *collideField, const int * const flagField, const t_procData * const procData, const int wallIdx){
//
// 	int endOuter = -1, endInner = -1, fixedValue = -1;
// 	p_setIterationParameters(&endOuter, &endInner, &fixedValue, procData, wallIdx);
// 	assert(endOuter != -1 && endInner != -1 && fixedValue != -1);
//
// 	int xlength[3] = {procData->xLength[0], procData->xLength[1], procData->xLength[2]};
//
// 	//variable needed at check whether it is an in-domain (FLUID) cell
// 	int maxValidIndex = Q*(xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2);
//
// 	//k - corresponds to the 'outer' value when computing the offset
// 	for(int k = 0; k <= endOuter; ++k){
// 		//j - corresponds to the 'inner' value
// 		for(int j = 0; j <= endInner; ++j){
//
// 			//flagField[xyz_offset] should only be a boundary cell (checked in p_setBounceBack)
// 			int xyz_offset = p_computeCellOffset(k, j, fixedValue, xlength, wallIdx);
// 			int current_cell_index = Q*xyz_offset;
// 			int boundaryType = flagField[xyz_offset];
//
// 			assert(boundaryType == NO_SLIP || boundaryType == MOVING_WALL);
//
// 			//PARALLEL boundaries are not treated when they occur (e.g. at shared edges)
// 			for(int i = 0; i < Q; ++i){
// 				int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
// 				int n_xyzoffset = p_computeNeighborCellOffset(k, j, fixedValue, c, xlength, wallIdx);
// 				int n_cell_index = Q*n_xyzoffset;
//
// 				//check (1) valid index: in case the direction of vector 'c' points to a non-existing cell
// 				//check (2) that the neighboring cell is a FLUID cell or PARALLEL_BOUNDARY (treated as inner)
// 				if(n_cell_index >= 0 && n_cell_index < maxValidIndex &&
// 						(flagField[n_xyzoffset] == FLUID || flagField[n_xyzoffset] == PARALLEL_BOUNDARY)
// 				){
// 					p_setBounceBack(collideField, procData->wallVelocity, boundaryType, i, current_cell_index, n_cell_index, c);
// 				}
// 			}
// 		}
// 	}
// }

#define VARIABLE -1 //symbol that is used in context of indices. It indicates which index is variable.
#define EXTRACT 0
#define INJECT 1
//TODO: (DL) possible define "9" for middle distribution - reduce magic numbers...

//TODO: (DL) Bug:
//As for the wall we take the entire wall for density. This is problematic, because if shared edges.
//At the moment the domain boundary edges should be treated correctly. However, the trouble comes with
//the subdomain shared edges...
//We need a special treatment for these (sub-domain) edges...

// |	*                            *   |
// |	*                            *   |
// |	*                          Z *   |
// |	******************************   |
// |_____________________________________|
//
// ______________________________________
// |O                                    |
// |   ******************************    |
// |   * X                          *    |
// |   *                            *    |
// |   *                            *    |
// |   *                            *    |

// Cell X needs the density of cell Z (though own streamfield boundary cell O).
// However, this communication is not yet done.

//The question is: Can we archive this with the already available SendRecvs or do need more?





//int extractInjectEdge(double buffer[], double * const collideField, t_iterParaEdge const * const iterPara,
		//t_procData const * const procData, const int index, const int injectFlag, const int densityFlag){

	//assert(injectFlag == EXTRACT || injectFlag == INJECT);

	//int endIndex, cellOffset;

	////Note: the corners of shared edges do not have to be treated
	////(in D3Q19 there are no distributions into the corners)
	//if(iterPara->x == VARIABLE){
		//endIndex = procData->xLength[0];
	//}else if(iterPara->y == VARIABLE){
		//endIndex = procData->xLength[1];
	//}else{
		//assert(iterPara->z == VARIABLE);
		//endIndex = procData->xLength[2];
	//}

	//int currentIndexBuff = 0;

	////z * (xlen*ylen) + y * (xlen) + x
	//for(int varIdx = 1; varIdx <= endIndex; ++varIdx){
		//if(iterPara->x == VARIABLE){
			//cellOffset = p_computeCellOffsetXYZ_Q(varIdx, iterPara->y, iterPara->z, procData->xLength);
		//}else if(iterPara->y == VARIABLE){
			//cellOffset = p_computeCellOffsetXYZ_Q(iterPara->x, varIdx, iterPara->z, procData->xLength);
		//}else{
			//assert(iterPara->z == VARIABLE);
			//cellOffset = p_computeCellOffsetXYZ_Q(iterPara->x, iterPara->y, varIdx, procData->xLength);
		//}

		//if(injectFlag){
			//if(densityFlag){
				//collideField[cellOffset+9] = buffer[currentIndexBuff++];
			//}else{
				//collideField[cellOffset+index] = buffer[currentIndexBuff++];
			//}

		//}else{
			//if(densityFlag){
				//c_computeNumDensity(&collideField[cellOffset], &buffer[currentIndexBuff]);
				//currentIndexBuff++;
			//}else{
				//buffer[currentIndexBuff++] = collideField[cellOffset+index];
			//}
		//}
	//}
	//return currentIndexBuff; //Will be used as bufferSize
//}

//void p_setBoundaryIterParameters(t_iterPara *const iterPara, t_procData const*const procData, const int direction){

	//iterPara->startOuter = 0;
	//iterPara->startInner = 0;

	//switch(direction/2){ //integer division to get the type of face (see enum in LBDefinitions.h)
	////---------------------------------------------
	////outer = Z, inner = X, Y fixed
	//case 0:
		//iterPara->endOuter   = procData->xLength[2]+1;
		//iterPara->endInner   = procData->xLength[0]+1;
		//iterPara->fixedValue = (direction == LEFT) ? 1 : procData->xLength[1];
		//break;

	////---------------------------------------------
	////outer = Y, inner = X, Z fixed
	//case 1:
		//iterPara->endOuter   = procData->xLength[1]+1;
		//iterPara->endInner   = procData->xLength[0]+1;
		//iterPara->fixedValue = (direction == BOTTOM) ? 1 : procData->xLength[2];
		//break;

	////---------------------------------------------
	////outer = Z, inner = Y, X fixed
	//case 2:
		//iterPara->endOuter   = procData->xLength[2]+1;
		//iterPara->endInner   = procData->xLength[1]+1;
		//iterPara->fixedValue = (direction == BACK) ? 1 : procData->xLength[0];
		//break;

	//default:
		//ERROR("Invalid direction occurred. This should never happen!");
	//}
//}

//void p_setEdgeIterParameters(t_iterParaEdge *const iterPara, t_procData const*const procData, const int edge, const int injectFlag){

	//assert(injectFlag == 0 || injectFlag == 1);

	////if inject: copy into boundary; else: copy from domain cells lying at edge
	//int start     = injectFlag ? 0 : 1; //if inject: start = 0, is at (ghost) boundary layer
	//int endOffset = injectFlag ? 1 : 0; //if inject: add endOffset +1 to be at (ghost) boundary layer

	//switch(edge/4){ //integer division - 12 edges, 3 cases
	//case 0: //edges 0-3 (horizontal, bottom)
		//iterPara->z = start;
		//if(edge == 0 || edge == 2){
			//iterPara->x = VARIABLE;
			//iterPara->y = (edge == 0) ? start : procData->xLength[1]+endOffset;
		//}else{ //edge == 1 || edge == 3
			//iterPara->y = VARIABLE;
			//iterPara->x = (edge == 1) ? procData->xLength[0]+endOffset : start;
		//}
		//break;

	//case 1: //edges 4-7 (vertical)
		//iterPara->z = VARIABLE;
		//iterPara->y = (edge == 4 || edge == 5) ? start : procData->xLength[1]+endOffset;
		//iterPara->x = (edge == 4 || edge == 7) ? start : procData->xLength[0]+endOffset;
		//break;

	//case 2: //edges 8-11 (horizontal, top)
		//iterPara->z = procData->xLength[2]+endOffset;
		//if(edge == 8 || edge == 10){
			//iterPara->x = VARIABLE;
			//iterPara->y = (edge == 8) ? start : procData->xLength[1]+endOffset;
		//}else{ //edge == 9 || edge == 11
			//iterPara->y = VARIABLE;
			//iterPara->x = (edge == 9) ? procData->xLength[0]+endOffset : start;
		//}
		//break;

	//default:
		//assert(edge>=0 && edge<=11);
	//}
	//assert(((iterPara->x < 0) + (iterPara->y < 0) + (iterPara->z < 0)) == 1); //only one is VARIABLE

	////TODO: (DL) delete when finalizing
	//// if(injectFlag == EXTRACT){
	//// 	printf("(EXTRACT) RANK: %i, x=%i, y=%i, z=%i, edge=%i \n", procData->rank, iterPara->x, iterPara->y, iterPara->z, edge);
	//// }else{
	//// 	printf("(INJECT) RANK: %i, x=%i, y=%i, z=%i, edge=%i \n", procData->rank, iterPara->x, iterPara->y, iterPara->z, edge);
	//// }
//}

//int p_assignSharedEdgeIndex(const int edge) {
	////Look for edges numbering in LBDefinitions.h
	//// edge 0: (0,-1,-1) 	-> [0]
	//// edge 1: (1,0,-1) 	-> [3]
	//// edge 2: (0,1,-1) 	-> [4]
	//// edge 3: (-1,0,-1) 	-> [1]

	//// edge 4: (-1,-1,0) 	-> [5]
	//// edge 5: (1,-1,0) 	-> [7]
	//// edge 6: (1,1,0) 		-> [13]
	//// edge 7: (-1,1,0) 	-> [11]

	//// edge 8: (0,-1,1) 	-> [14]
	//// edge 9: (1,0,1) 		-> [17]
	//// edge 10: (0,1,1) 	-> [18]
	//// edge 11: (-1,0,1) 	-> [15]
	//static int indices[12] = {0,3,4,1,
							  //5,7,13,11,
							  //14,17,18,15};
	//return indices[edge];
//}

//void treatPeriodicWall(int const * const flagField, double *const collideField, double *const sendBuffer, double *const readBuffer,
	//const t_procData * const procData, const int procWall, const int opponentWall, const int densityFlag){

	//t_iterPara  iterPara;
	//int 		indexIn[5], indexOut[5], bufferSize, commRank;

	//p_assignIndices(procWall,     indexOut);
	//p_assignIndices(opponentWall, indexIn);

	//p_setBoundaryIterParameters(&iterPara, procData, procWall);

	////always excluding shared edges (treat differently&separately)
	//bufferSize = densityFlag ? (iterPara.endOuter+1) * (iterPara.endInner+1) :  5 * (iterPara.endOuter+1) * (iterPara.endInner+1);

	//// printf("startInner=%i, startOuter=%i, endInner=%i, endOuter=%i \n",
	//// iterPara.startInner, iterPara.startOuter, iterPara.endInner, iterPara.endOuter);
	//if(densityFlag){
		//extractDensity(flagField, sendBuffer, collideField, &iterPara, procData, procWall);
	//}else{
		//extract(sendBuffer, collideField, &iterPara, procData, procWall, indexOut);
	//}

	//commRank = procData->periodicNeighbours[procWall];
	//assert(commRank >= 0);

	//// if(procData->rank == 1 && procWall == TOP){
	//// 	for(int i = 0; i < bufferSize; i++){
	//// 		printf("SEND: %f \n", sendBuffer[i]);
	//// 	}
	//// }

	////tag is chosen to be 1 that it is different from parallel_boundary (0) and edges (2)
	//MPI_Sendrecv(sendBuffer, bufferSize, MPI_DOUBLE, commRank, 1, readBuffer,
		//bufferSize, MPI_DOUBLE, commRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	//// if(procData->rank == 0 && procWall == BOTTOM){
	//// 	for(int i = 0; i < bufferSize; i++){
	//// 		printf("READ: %f \n", readBuffer[i]);
	//// 	}
	//// }

	//if(densityFlag){
		//const int indexIn = 9;
		//inject(readBuffer, collideField, &iterPara, procData, procWall, &indexIn, 1);
	//}else{
		//inject(readBuffer, collideField, &iterPara, procData, procWall, indexIn, nrDistSwap);
	//}
//}

//void treatPeriodicEdge(double *collideField, double *const sendBuffer, double *const readBuffer,
	//const t_procData * const procData, const int procEdge, const int opponentEdge, const int densityFlag){

	//t_iterParaEdge iterParaEdge;
	////Note: in edge case only one index has to be copied
	//int indexOutEdge, indexInEdge, bufferSize, commRank;

	//indexOutEdge = densityFlag ? -1 : p_assignSharedEdgeIndex(procEdge);
	//indexInEdge =  densityFlag ? -1 : p_assignSharedEdgeIndex(opponentEdge);

	////printf("RANK: %i: indexOutEdge = %i, indexInEdge = %i \n", procData->rank, procEdge, opponentEdge);
	//p_setEdgeIterParameters(&iterParaEdge, procData, procEdge, EXTRACT);

	//// if(procData->rank == 0){
	//// 	printf("(EXTRACT) procData->rank=%i, procEdge=%i, opponentEdge=%i ",procData->rank, procEdge, opponentEdge);
	//// 	printf("iterParaEdge.x=%i, iterParaEdge.y=%i, iterParaEdge.z=%i \n", iterParaEdge.x, iterParaEdge.y, iterParaEdge.z);
	//// }

	////using buffers from parallel boundaries
	//bufferSize = extractInjectEdge(sendBuffer, collideField, &iterParaEdge, procData, indexOutEdge, EXTRACT, densityFlag);

	//commRank = procData->periodicEdgeNeighbours[procEdge]; //communication rank
	//assert(commRank >= 0);

	////tag is chosen to be 2 that it is different from parallel_boundary (0) and periodic_boundary (1)
	//MPI_Sendrecv(sendBuffer, bufferSize, MPI_DOUBLE, commRank, 2, readBuffer,
			//bufferSize, MPI_DOUBLE, commRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	//p_setEdgeIterParameters(&iterParaEdge, procData, procEdge, INJECT);

	//// if(procData->rank == 0){
	//// 	printf("(INJECT) procData->rank=%i, procEdge=%i, opponentEdge=%i ",procData->rank, procEdge, opponentEdge);
	//// 	printf("iterParaEdge.x=%i, iterParaEdge.y=%i, iterParaEdge.z=%i \n", iterParaEdge.x, iterParaEdge.y, iterParaEdge.z);
	//// }
	//extractInjectEdge(readBuffer, collideField, &iterParaEdge, procData, indexInEdge, INJECT, densityFlag);
//}

//void treatPeriodicWallNoComm(double *collideField, const int wall1, const int wall2, t_procData const*const procData, const int densityFlag){
	//t_iterPara  iterPara1 = {0}, iterPara2 = {0}; //default initialization to avoid "may be uninitialized" warning (=error)

	//int index1[5], index2[5];
	//p_setBoundaryIterParameters(&iterPara1, procData, wall1);
	//p_setBoundaryIterParameters(&iterPara2, procData, wall2);

	//if(! densityFlag){
		//p_assignIndices(wall1, index1);
		//p_assignIndices(wall2, index2);
	//}

	////only fixed value should be different
	//assert(iterPara1.startInner == iterPara2.startInner && iterPara1.startOuter == iterPara2.startOuter &&
		   //iterPara1.endInner == iterPara2.endInner && iterPara1.endOuter == iterPara2.endOuter);

	//int currentIndexFieldIn1 = -1, currentIndexFieldOut1 = -1,
		//currentIndexFieldIn2 = -1, currentIndexFieldOut2 = -1;//initially invalid assignment

	//static int shiftFixedValue[6] = {-1, 1, 1, -1, 1, -1}; //static to allocate/assign only once

	////k - corresponds to the 'outer' value when computing the offset
	//for(int k = iterPara1.startOuter; k <= iterPara1.endOuter; ++k){
		////j - corresponds to the 'inner' value
		//for(int j = iterPara1.startInner; j <= iterPara1.endInner; ++j){

			//currentIndexFieldIn1 = Q*p_computeCellOffset(k, j, iterPara1.fixedValue+shiftFixedValue[wall1], procData->xLength, wall1);
			//currentIndexFieldOut1 = Q*p_computeCellOffset(k, j, iterPara1.fixedValue, procData->xLength, wall1);

			//currentIndexFieldIn2 = Q*p_computeCellOffset(k, j, iterPara2.fixedValue+shiftFixedValue[wall2], procData->xLength, wall2);
			//currentIndexFieldOut2 = Q*p_computeCellOffset(k, j, iterPara2.fixedValue, procData->xLength, wall2);

			//// Out of bounds check (only 'out' values at the moment)
			//assert(currentIndexFieldOut1 < Q*(procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2)
			//&& currentIndexFieldOut1 >= 0);
			//assert(currentIndexFieldOut2 < Q*(procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2)
			//&& currentIndexFieldOut2 >= 0);

			//if(densityFlag){
				////write the density to "In" form the cell "Out"
				//c_computeNumDensity(&collideField[currentIndexFieldOut2], &collideField[currentIndexFieldIn1+9]);
				//c_computeNumDensity(&collideField[currentIndexFieldOut1], &collideField[currentIndexFieldIn2+9]);
			//}else{
				////Do copy operations:
				//for (int i = 0; i < nrDistSwap; i++) {
					//collideField[currentIndexFieldIn1+index1[i]] = collideField[currentIndexFieldOut2+index1[i]];
					//collideField[currentIndexFieldIn2+index2[i]] = collideField[currentIndexFieldOut1+index2[i]];
				//}
			//}
		//}
	//}
//}

//void treatPeriodicEdgeNoComm(double *collideField, const t_procData * const procData, const int edge1, const int edge2, const int densityFlag){

	//int endIndex;
	//t_iterParaEdge iterParaEdgeIn1, iterParaEdgeOut1, iterParaEdgeIn2, iterParaEdgeOut2;
	//int index1, index2;

	//index1 = p_assignSharedEdgeIndex(edge1);
	//index2 = p_assignSharedEdgeIndex(edge2);
	//assert(index1 == Q-index2-1);

	////fixed values change, iterParameters have to have the same value set to VARIABLE
	//p_setEdgeIterParameters(&iterParaEdgeIn1, procData, edge1, INJECT);
	//p_setEdgeIterParameters(&iterParaEdgeOut1, procData, edge1, EXTRACT);

	//p_setEdgeIterParameters(&iterParaEdgeIn2, procData, edge2, INJECT);
	//p_setEdgeIterParameters(&iterParaEdgeOut2, procData, edge2, EXTRACT);

	////Note: the corners of shared edges do not have to be treated
	////(in D3Q19 there are no distributions into the corners)
	//if(iterParaEdgeIn1.x == VARIABLE){
		//assert(iterParaEdgeOut1.x == VARIABLE && iterParaEdgeIn2.x == VARIABLE && iterParaEdgeOut2.x == VARIABLE);
		//endIndex = procData->xLength[0];
	//}else if(iterParaEdgeIn1.y == VARIABLE){
		//assert(iterParaEdgeOut1.y == VARIABLE && iterParaEdgeIn2.y == VARIABLE && iterParaEdgeOut2.y == VARIABLE);
		//endIndex = procData->xLength[1];
	//}else{
		//assert(iterParaEdgeIn1.z == VARIABLE && iterParaEdgeOut1.z == VARIABLE && iterParaEdgeIn2.z == VARIABLE && iterParaEdgeOut2.z == VARIABLE);
		//endIndex = procData->xLength[2];
	//}

	//int cellOffsetIn1, cellOffsetOut1, cellOffsetIn2, cellOffsetOut2;
	////z * (xlen*ylen) + y * (xlen) + x
	//for(int varIdx = 1; varIdx <= endIndex; ++varIdx){
		//if(iterParaEdgeIn1.x == VARIABLE){
			//cellOffsetIn1 	= p_computeCellOffsetXYZ_Q(varIdx, iterParaEdgeIn1.y, iterParaEdgeIn1.z, procData->xLength);
			//cellOffsetOut1	= p_computeCellOffsetXYZ_Q(varIdx, iterParaEdgeOut1.y, iterParaEdgeOut1.z, procData->xLength);
			//cellOffsetIn2 	= p_computeCellOffsetXYZ_Q(varIdx, iterParaEdgeIn2.y, iterParaEdgeIn2.z, procData->xLength);
			//cellOffsetOut2 	= p_computeCellOffsetXYZ_Q(varIdx, iterParaEdgeOut2.y, iterParaEdgeOut2.z, procData->xLength);

		//}else if(iterParaEdgeIn1.y == VARIABLE){
			//cellOffsetIn1 	= p_computeCellOffsetXYZ_Q(iterParaEdgeIn1.x, varIdx, iterParaEdgeIn1.z, procData->xLength);
			//cellOffsetOut1	= p_computeCellOffsetXYZ_Q(iterParaEdgeOut1.x, varIdx, iterParaEdgeOut1.z, procData->xLength);
			//cellOffsetIn2 	= p_computeCellOffsetXYZ_Q(iterParaEdgeIn2.x, varIdx, iterParaEdgeIn2.z, procData->xLength);
			//cellOffsetOut2 	= p_computeCellOffsetXYZ_Q(iterParaEdgeOut2.x, varIdx, iterParaEdgeOut2.z, procData->xLength);
		//}else{
			//cellOffsetIn1 	= p_computeCellOffsetXYZ_Q(iterParaEdgeIn1.x, iterParaEdgeIn1.y, varIdx, procData->xLength);
			//cellOffsetOut1	= p_computeCellOffsetXYZ_Q(iterParaEdgeOut1.x, iterParaEdgeOut1.y, varIdx, procData->xLength);
			//cellOffsetIn2 	= p_computeCellOffsetXYZ_Q(iterParaEdgeIn2.x, iterParaEdgeIn2.y, varIdx, procData->xLength);
			//cellOffsetOut2 	= p_computeCellOffsetXYZ_Q(iterParaEdgeOut2.x, iterParaEdgeOut2.y, varIdx, procData->xLength);
		//}

		//if(densityFlag){
			////+9 is the distribution at the center (where we save the number density)
			//c_computeNumDensity(&collideField[cellOffsetOut2], &collideField[cellOffsetIn1+9]);
			//c_computeNumDensity(&collideField[cellOffsetOut1], &collideField[cellOffsetIn2+9]);
		//}else{
			//collideField[cellOffsetIn1+index1] = collideField[cellOffsetOut2+index1];
			//collideField[cellOffsetIn2+index2] = collideField[cellOffsetOut1+index2];
		//}
	//}
//}

//void treatComponentBoundary(t_component *c, int const*const flagField, t_procData const*const procData, double **sendBuffer, double **readBuffer, const int densityFlag){
	//for (int i = 0; i < numComp; ++i) {
        //treatBoundary(flagField, c[i].collideField, procData, sendBuffer, readBuffer, densityFlag);
    //}
//}

////TODO: (TKS) Remove flagField if not in use when finished
//void treatBoundary(int const*const flagField, double *const collideField, const t_procData * const procData, double **sendBuffer, double **readBuffer, const int densityFlag){
	//assert(densityFlag == 0 || densityFlag == 1);

	//// Handle of "real" boundaries:
	//// const int NO_NEIGHBOUR = -2; // equals the MPI_PROC_NULL = -2
	//// for(int wall=LEFT; wall<=BACK; ++wall){ //see LBDefinitions for order of walls
	//// 	if(procData->neighbours[wall] == NO_NEIGHBOUR){
	//// 		// printWallEnum(wall);
	//// 		p_treatSingleWall(collideField, flagField, procData, wall);
	//// 	}
	//// }

	//for(int wall=LEFT; wall<=BACK; wall+=2){ //see LBDefinitions for order of walls

		//int firstNeighbour = procData->periodicNeighbours[wall];
		//int secondNeighbour = procData->periodicNeighbours[wall+1];

		//if(firstNeighbour == procData->rank && secondNeighbour == procData->rank){
			////Case when there is only one proc, then no communication is required
			//treatPeriodicWallNoComm(collideField, wall, wall+1, procData, densityFlag);
		//}
		//else{
			//if(firstNeighbour != MPI_PROC_NULL){
				//treatPeriodicWall(flagField,collideField, sendBuffer[wall], readBuffer[wall],
					//procData, wall, wall+1, densityFlag);
			//}

			//if(secondNeighbour != MPI_PROC_NULL){
				//treatPeriodicWall(flagField,collideField, sendBuffer[wall+1], readBuffer[wall+1],
					//procData, wall+1, wall, densityFlag);
			//}
		//}
	//}

	//static const int edge1[6] = {0,1,2,3,4,5};
	//static const int edge2[6] = {10,11,8,9,6,7}; //opposite edge to edge1; see numbering in LBDefinitions.h
	////TODO: (DL) this is not nice, think about alternatives. For now there is guaranteed enough space.
	//double *validSendBuffer, *validReadBuffer;

	//for(int idx = 0; idx < 6; ++idx){

		//int firstNeighbour = procData->periodicEdgeNeighbours[edge1[idx]];
		//int secondNeighbour = procData->periodicEdgeNeighbours[edge2[idx]];

		////TODO: (DL) make a better solution for buffer business
		//if(edge1[idx] == 0 || edge1[idx] == 2 || edge1[idx] == 4 || edge1[idx] == 5){
			//validSendBuffer = sendBuffer[LEFT];
			//validReadBuffer = readBuffer[LEFT];
		//}else{ //(edge1[idx] == 1 || edge1[idx] == 3){
			//assert(edge1[idx] == 1 || edge1[idx] == 3);
			//validSendBuffer = sendBuffer[TOP];
			//validReadBuffer = readBuffer[TOP];
		//}

		//if(firstNeighbour == procData->rank && secondNeighbour == procData->rank){
			////Case when there is only one proc, then no communication is required
			//treatPeriodicEdgeNoComm(collideField, procData, edge1[idx], edge2[idx], densityFlag);
		//}else{
			//if(firstNeighbour != MPI_PROC_NULL){
				//treatPeriodicEdge(collideField, validSendBuffer, validReadBuffer,
					//procData, edge1[idx], edge2[idx], densityFlag);
			//}

			//if(secondNeighbour != MPI_PROC_NULL){
				//treatPeriodicEdge(collideField, validSendBuffer, validReadBuffer,
					//procData, edge2[idx], edge1[idx], densityFlag);
			//}
		//}
	//}
//}

void treatBoundary(t_component *c, int xlength){
	for (int i = 0; i < NUMCOMP; i++) {
		treatWallPeriodic(&c[i], LEFT, xlength);
		treatWallPeriodic(&c[i], RIGHT, xlength);
		treatWallPeriodic(&c[i], TOP, xlength);
		treatWallPeriodic(&c[i], BOTTOM, xlength);
		treatWallPeriodic(&c[i], FRONT, xlength);
		treatWallPeriodic(&c[i], BACK, xlength);
	}
}

// This will be useful for computing otherSideIdx
// nbhR = (myRank+1)%numRanks;
// nbhL = (myRank-1+numRanks)%numRanks;
void treatWallPeriodic(t_component * c, int direction, int xlength) {
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
			int idx = computeCellOffset(k, j, fixedValueIdx, direction, xlength);
			int idx_Q = Q*idx;
			int otherSideIdx = computeCellOffset(k, j, fixedValueOtherIdx, direction, xlength);
			int otherSideIdx_Q = Q*otherSideIdx;

			c->rho[otherSideIdx] = c->rho[idx];
			for (int i = 0; i < Q; i++) {
                //TODO: (TKS) Commenting out stream and collide gives succesfull test.
				//c->collideField[otherSideIdx_Q+i] = c->collideField[idx_Q+i];
				//c->streamField[otherSideIdx_Q+i] = c->streamField[idx_Q+i];
				c->feq[otherSideIdx_Q+i] = c->feq[idx_Q+i];
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
