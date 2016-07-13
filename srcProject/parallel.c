#include "parallel.h"
#include "LBDefinitions.h"

void initialiseMPI(int *rank, int *numRanks, int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,numRanks);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);
}

void broadcastValues(int rank, int *xlength, double* rhoFluct, t_component *c, double G[numComp][numComp],
	int *procsPerAxis, int *timesteps, int *timestepsPerPlotting, int* timestepDouble, int* timestepMax) {

	MPI_Bcast(xlength, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(timestepsPerPlotting, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(timestepDouble, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(timestepMax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(procsPerAxis, 3, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(G, numComp*numComp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(rhoFluct, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int i = 0; i < numComp; i++) {
		MPI_Bcast(&c[i].tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&c[i].m, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&c[i].rhoRef, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	 	MPI_Bcast(&c[i].psiFctCode, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

}

/*
* Divides the entire domain into sub-domains. The function sets the processor specific data to determine the size
* of the sub-domain and determines the neighbours.
*/
void domainDecompositionAndNeighbors(t_procData *const procData, const int * const xlength, int const * const procsPerAxis) {

	#define M_NBGH procData->neighbours

	int procPos[3] = {0,0,0};

	// Get the position of the process based on its rank
	p_rankToPos(procsPerAxis,procData->rank,procPos);

	// Compute the subdomain size and save it into procData
	int baseLength[3] = {xlength[0]/procsPerAxis[0], xlength[1]/procsPerAxis[1], xlength[2]/procsPerAxis[2]};

	// If the proc is at the end of some axis then add the remaining length
	procData->xLength[0] =  (procPos[0] == procsPerAxis[0]-1)?
                            xlength[0]-(procsPerAxis[0]-1)*baseLength[0]:baseLength[0];

	procData->xLength[1] = (procPos[1] == procsPerAxis[1]-1)?
                           xlength[1]-(procsPerAxis[1]-1)*baseLength[1]:baseLength[1];

	procData->xLength[2] = (procPos[2] == procsPerAxis[2]-1)?
                           xlength[2]-(procsPerAxis[2]-1)*baseLength[2]:baseLength[2];

	// Set neighbours for parallel boundaries
	// Decide whether it is a ghost boundary (= MPI_PROC_NULL) or a parallel boundary (assign rank of neighbor)
    procData->neighbours[LEFT]	 = procPos[0]
                                 + (procPos[1]-1+procsPerAxis[1])%procsPerAxis[1]*procsPerAxis[0]
                                 + procPos[2]*procsPerAxis[0]*procsPerAxis[1];

    procData->neighbours[RIGHT]	 = procPos[0]
                                 + (procPos[1]+1)%procsPerAxis[1]*procsPerAxis[0]
                                 + procPos[2]*procsPerAxis[0]*procsPerAxis[1];

    procData->neighbours[TOP]	 = procPos[0]
                                 + procPos[1]*procsPerAxis[0]
                                 + (procPos[2]+1)%procsPerAxis[2]*procsPerAxis[0]*procsPerAxis[1];

    procData->neighbours[BOTTOM] = procPos[0]
                                 + procPos[1]*procsPerAxis[0]
                                 + (procPos[2]-1+procsPerAxis[2])%procsPerAxis[2]*procsPerAxis[0]*procsPerAxis[1];

    procData->neighbours[FRONT]	 = (procPos[0]+1)%procsPerAxis[0]
                                 + procPos[1]*procsPerAxis[0]
                                 + procPos[2]*procsPerAxis[0]*procsPerAxis[1];

    procData->neighbours[BACK]	 = (procPos[0]-1+procsPerAxis[0])%procsPerAxis[0]
                                 + procPos[1]*procsPerAxis[0]
                                 + procPos[2]*procsPerAxis[0]*procsPerAxis[1];


	// If self neighbor then set neighbor to MPI_PROC_NULL
	for (int i = 0; i < 6; i++) {
		if (procData->neighbours[i] == procData->rank) {
			procData->neighbours[i] = MPI_PROC_NULL;
		}
	}
}

/*
* Initialize sendBuffer and readBuffer. They are only initialized in a certain
* direction when there is a valid neighbor, i.e. no domain boundary
*/
void initialiseBuffers(double *sendBuffer[6], double *readBuffer[6], int const * const xlength,
	int const * const neighbours, int procBufferSize[3]) {

    const int db = sizeof(double); //double in bytes


    //---------------------------------------------
    //      NOTE
    //---------------------------------------------
    // In the current implementation all buffers
    // are initialized, if we had other BC than
    // periodic, e.g. noslip,  this should not be
    // done to save memory.
    //---------------------------------------------

    // XZ inner domain (no edges included)
    int bufferSize		= nrDistSwap*(xlength[0]*(xlength[2]));
    procBufferSize[0] 	= bufferSize; //Valid for left and right

	sendBuffer[LEFT] 	= (double *) malloc(bufferSize*db);
	sendBuffer[RIGHT] 	= (double *) malloc(bufferSize*db);
	readBuffer[LEFT] 	= (double *) malloc(bufferSize*db);
	readBuffer[RIGHT] 	= (double *) malloc(bufferSize*db);


    // XY plane including edges at the boundary to left/right
    bufferSize 			= nrDistSwap*(xlength[0]*(xlength[1]+2));
    procBufferSize[1] 	= bufferSize; //Valid for top and bottom

    sendBuffer[TOP] 	= (double *) malloc(bufferSize*db);
    sendBuffer[BOTTOM] 	= (double *) malloc(bufferSize*db);
    readBuffer[TOP] 	= (double *) malloc(bufferSize*db);
    readBuffer[BOTTOM] 	= (double *) malloc(bufferSize*db);


    // YZ plane including all edges
    bufferSize 			= nrDistSwap*((xlength[1]+2)*(xlength[2]+2));
    procBufferSize[2] 	= bufferSize; //Valid for front and back

	sendBuffer[FRONT] 	= (double *) malloc(bufferSize*db);
	sendBuffer[BACK] 	= (double *) malloc(bufferSize*db);
	readBuffer[FRONT] 	= (double *) malloc(bufferSize*db);
	readBuffer[BACK] 	= (double *) malloc(bufferSize*db);
}

/*
* Wrapper around communicate to communicate each component
*/
void communicateComponents(double** sendBuffer, double**readBuffer, t_component *c, t_procData const * const procData){
    for (int i = 0; i < numComp; ++i) {
        communicate(sendBuffer, readBuffer, c[i].field, procData, 0);
		communicate(sendBuffer, readBuffer, c[i].feq, procData, 1);
		communicateDensity(sendBuffer, readBuffer, c[i].rho, procData);
    }
}


/*
 * Main algorithm to carry out the communication between the processes.
 */
void communicate(double** sendBuffer, double**readBuffer, double* field_new, t_procData const * const procData, int tag){
    //Run extract, swap, inject for all sides and cells.
    int 		index1[5], index2[5];
    t_iterPara  iterPara1, iterPara2;

	// Iterate through directions (left/right, top/bottom, front/back)
    for (int direction = LEFT; direction <= BACK; direction+=2) {

        const int face1 = direction;   //LEFT,  TOP,    FRONT
        const int face2 = direction+1; //RIGHT, BOTTOM, BACK

		// Set iteration parameters for both faces of a direction
        p_setCommIterationParameters(&iterPara1, procData, face1);
        p_setCommIterationParameters(&iterPara2, procData, face2);

		// Assign the indices for the 5 distributions that go out of the two faces of the current direction
        p_assignIndices(face1, index1);
        p_assignIndices(face2, index2);

        if(procData->neighbours[face1] != MPI_PROC_NULL && procData->neighbours[face2] != MPI_PROC_NULL) {
		    // Extract the distributions from field_new to the send buffer
			extract(sendBuffer[face1], field_new, &iterPara1, procData, face1, index1);
			extract(sendBuffer[face2], field_new, &iterPara2, procData, face2, index2);

            // Swap the distributions
            swap(sendBuffer, readBuffer, procData, direction, tag);

            // Inject the distributions from read buffer to collide field
            //NOTE: opposite index get injected from the readBuffer
            inject(readBuffer[face1], field_new, &iterPara1, procData, face1, index2);
            inject(readBuffer[face2], field_new, &iterPara2, procData, face2, index1);
		}
        else if (procData->neighbours[face1] == MPI_PROC_NULL && procData->neighbours[face2] == MPI_PROC_NULL) {
			swapNoComm(field_new, &iterPara1, &iterPara2, procData, direction, index1, index2);
		}
        else {
            ERROR("No communication of distributions happened!");
        }
    }
}

/*
 * Copy distributions into sendBuffer.
 */
void extract(double sendBuffer[], double const * const field_new, t_iterPara const * const iterPara, t_procData const * const procData,
             const int direction, int const * const index){

    int currentIndexField = -1; //initially invalid assignment
    int currentIndexBuff  = 0;  //simple counter

    //k - corresponds to the 'outer' value when computing the offset
    for(int k = iterPara->startOuter; k <= iterPara->endOuter; ++k){
        //j - corresponds to the 'inner' value
        for(int j = iterPara->startInner; j <= iterPara->endInner; ++j){

            currentIndexField  = Q*p_computeCellOffset(k, j, iterPara->fixedValue, procData->xLength, direction);

			// Out of bounds check
            assert(currentIndexField < Q*(procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2)
                    && currentIndexField >= 0);

            for (int i = 0; i < 5; i++) {
				// Out of bounds check
				assert(currentIndexBuff < procData->bufferSize[direction/2] && currentIndexBuff >=0);
                sendBuffer[currentIndexBuff++] = field_new[currentIndexField+index[i]];
            }
        }
    }
}


/*
 * Swap distributions with neighboring processes.
 */
void swap(double** sendBuffer, double** readBuffer, t_procData const * const procData, const int direction, int tag){

    //temp variables for readabillity.
	const int bufferSize = procData->bufferSize[direction/2];
    const int proc1 = procData->neighbours[direction];
    const int proc2 = procData->neighbours[direction+1];

    //MPI Sendrecv API
	//int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //                 int dest, int sendtag, void *recvbuf, int recvcount,
    //                 MPI_Datatype recvtype, int source, int recvtag,
    //                 MPI_Comm comm, MPI_Status *status)
    //Send proc1 receive proc1

    if(procData->neighbours[direction] < procData->rank){

        MPI_Sendrecv(sendBuffer[direction], bufferSize , MPI_DOUBLE, proc1, tag, readBuffer[direction],
            bufferSize, MPI_DOUBLE, proc1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendBuffer[direction+1], bufferSize , MPI_DOUBLE, proc2, tag, readBuffer[direction+1],
            bufferSize, MPI_DOUBLE, proc2, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Sendrecv(sendBuffer[direction+1], bufferSize , MPI_DOUBLE, proc2, tag, readBuffer[direction+1],
            bufferSize, MPI_DOUBLE, proc2, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendBuffer[direction], bufferSize , MPI_DOUBLE, proc1, tag, readBuffer[direction],
            bufferSize, MPI_DOUBLE, proc1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

/*
 * Copy read buffer into ghost layer of field_new.
 */
void inject(double const * const readBuffer, double* field_new, t_iterPara const * const iterPara, t_procData const * const procData,
            const int direction, int const * const index){

    int currentIndexField = -1; //Initially assign to invalid
    int currentIndexBuff = 0;   //Simple counter


    static int shiftFixedValue[6] = {-1, 1, 1, -1, 1, -1}; //static to allocate/assign only once

    //change value also in outer scope okay, because not needed afterwards
    int newFixedValue = iterPara->fixedValue + shiftFixedValue[direction];

    //k - corresponds to the 'outer' value when computing the offset
    for(int k = iterPara->startOuter; k <= iterPara->endOuter; ++k){
        //j - corresponds to the 'inner' value
        for(int j = iterPara->startInner; j <= iterPara->endInner; ++j){

            currentIndexField  = Q*p_computeCellOffset(k, j, newFixedValue, procData->xLength, direction);

            //out of bound checks
            assert(currentIndexField < Q*(procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2)
                    && currentIndexField >= 0);

            for (int i = 0; i < 5; i++) {
				//out of bound checks
				assert(currentIndexBuff < procData->bufferSize[direction/2] && currentIndexBuff >=0);
                field_new[currentIndexField+index[i]] = readBuffer[currentIndexBuff++];
            }
        }
    }
}

/*For periodic boundaries where one process can be its own neihbour. Does the same as swap from above*/
void swapNoComm(double* field_new, t_iterPara const * const iterPara1,t_iterPara const * const iterPara2,
	t_procData const * const procData, const int direction, int const * const index1, int const * const index2) {

	int fromFixedValue1, toFixedValue1, fromFixedValue2, toFixedValue2;
	if (iterPara1->fixedValue == 1) {
		fromFixedValue1 = iterPara1->fixedValue;
		toFixedValue1   = iterPara2->fixedValue+1;
		fromFixedValue2 = iterPara2->fixedValue;
		toFixedValue2   = iterPara1->fixedValue-1;
	} else {
		fromFixedValue1 = iterPara1->fixedValue;
		toFixedValue1   = iterPara2->fixedValue-1;
		fromFixedValue2 = iterPara2->fixedValue;
		toFixedValue2   = iterPara1->fixedValue+1;
	}
	for (int k = iterPara1->startOuter; k <= iterPara1->endOuter; k++) {
		for (int j = iterPara1->startInner; j <= iterPara1->endInner; j++) {

			int fromIdx1 = Q*p_computeCellOffset(k, j, fromFixedValue1, procData->xLength, direction);
			int toIdx1   = Q*p_computeCellOffset(k, j, toFixedValue1, procData->xLength, direction);
			int fromIdx2 = Q*p_computeCellOffset(k, j, fromFixedValue2, procData->xLength, direction);
			int toIdx2   = Q*p_computeCellOffset(k, j, toFixedValue2, procData->xLength, direction);

			for (int i = 0; i < 5; i++) {
				field_new[toIdx1+index1[i]] = field_new[fromIdx1+index1[i]];
				field_new[toIdx2+index2[i]] = field_new[fromIdx2+index2[i]];
			}
		}
	}
}

/*Communicates densities to neighbouring processes*/
void communicateDensity(double** sendBuffer, double**readBuffer, double* rho,
	t_procData const * const procData) {

    t_iterPara  iterPara1, iterPara2;

	// Iterate through directions (left/right, top/bottom, front/back)
    for (int direction = LEFT; direction <= BACK; direction+=2) {

        const int face1 = direction;   //LEFT,  TOP,    FRONT
        const int face2 = direction+1; //RIGHT, BOTTOM, BACK

		// Set iteration parameters for both faces of a direction
        p_setCommIterationParameters(&iterPara1, procData, face1);
        p_setCommIterationParameters(&iterPara2, procData, face2);

		int bufferSize1 = 0, bufferSize2 = 0;

		// Extract the distributions from field_new to the send buffer
        if(procData->neighbours[face1] != MPI_PROC_NULL && procData->neighbours[face2] != MPI_PROC_NULL) {
			bufferSize1 = extractDensity(sendBuffer[face1], rho, &iterPara1, procData, face1);
			bufferSize2 = extractDensity(sendBuffer[face2], rho, &iterPara2, procData, face2);

			assert(bufferSize1 == procData->bufferSize[direction/2]/nrDistSwap);
			assert(bufferSize2 == procData->bufferSize[direction/2]/nrDistSwap);

            // Swap the distributions
            swapDensity(sendBuffer, readBuffer, procData, direction, bufferSize1, bufferSize2);

            // Inject the distributions from read buffer to collide field
            //NOTE: opposite index get injected from the readBuffer
            injectDensity(readBuffer[face1], rho, &iterPara1, procData, face1);
            injectDensity(readBuffer[face2], rho, &iterPara2, procData, face2);
		}
        else if (procData->neighbours[face1] == MPI_PROC_NULL && procData->neighbours[face2] == MPI_PROC_NULL) {
			swapNoCommDensity(rho, &iterPara1, &iterPara2, procData, direction);
		}
        else{
            ERROR("No communication of densities happened!");

        }
    }
}


int extractDensity(double sendBuffer[], double const * const rho,
	t_iterPara const * const iterPara, t_procData const * const procData, const int direction) {

	int currentCellIndex = -1; //Initially assign to invalid
    int currentIndexBuff  = 0; //simple counter

    //k - corresponds to the 'outer' value when computing the offset
    for(int k = iterPara->startOuter; k <= iterPara->endOuter; ++k){
        //j - corresponds to the 'inner' value
        for(int j = iterPara->startInner; j <= iterPara->endInner; ++j){

			currentCellIndex = p_computeCellOffset(k, j, iterPara->fixedValue, procData->xLength, direction);

			// Out of bounds check
            assert(currentCellIndex < (procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2)
                    && currentCellIndex >= 0);

			sendBuffer[currentIndexBuff++] = rho[currentCellIndex];
        }
    }
	return currentIndexBuff;
}


void swapDensity(double** sendBuffer, double** readBuffer, t_procData const * const procData,
	const int direction, const int bufferSize1, const int bufferSize2){

    //temp variables for readabillity.
    const int proc1 = procData->neighbours[direction];
    const int proc2 = procData->neighbours[direction+1];

    //MPI Sendrecv API
	//int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //                 int dest, int sendtag, void *recvbuf, int recvcount,
    //                 MPI_Datatype recvtype, int source, int recvtag,
    //                 MPI_Comm comm, MPI_Status *status)
    //Send proc1 receive proc1
    assert(bufferSize1 != 0 && bufferSize2 != 0);
    if(procData->neighbours[direction] < procData->rank){
        MPI_Sendrecv(sendBuffer[direction], bufferSize1 , MPI_DOUBLE, proc1, 0, readBuffer[direction],
            bufferSize1, MPI_DOUBLE, proc1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(sendBuffer[direction+1], bufferSize2 , MPI_DOUBLE, proc2, 0, readBuffer[direction+1],
            bufferSize2, MPI_DOUBLE, proc2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Sendrecv(sendBuffer[direction+1], bufferSize2 , MPI_DOUBLE, proc2, 0, readBuffer[direction+1],
            bufferSize2, MPI_DOUBLE, proc2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(sendBuffer[direction], bufferSize1, MPI_DOUBLE, proc1, 0, readBuffer[direction],
            bufferSize1, MPI_DOUBLE, proc1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

int injectDensity(double const * const readBuffer, double* rho,
	t_iterPara const * const iterPara, t_procData const * const procData, const int direction){

	int currentFlagIndex = -1;
    int currentIndexBuff = 0;

    static int shiftFixedValue[6] = {-1, 1, 1, -1, 1, -1}; //static to allocate/assign only once

    //change value also in outer scope okay, because not needed afterwards
    int newFixedValue = iterPara->fixedValue + shiftFixedValue[direction];

    //k - corresponds to the 'outer' value when computing the offset
    for(int k = iterPara->startOuter; k <= iterPara->endOuter; ++k){
        //j - corresponds to the 'inner' value
        for(int j = iterPara->startInner; j <= iterPara->endInner; ++j){
			currentFlagIndex = p_computeCellOffset(k, j, newFixedValue, procData->xLength, direction);

            //out of bound checks
            assert(currentFlagIndex < (procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2)
                    && currentFlagIndex >= 0);

			assert(currentIndexBuff < procData->bufferSize[direction/2] && currentIndexBuff >=0);

            rho[currentFlagIndex] = readBuffer[currentIndexBuff++];
        }
    }
	return currentIndexBuff;
}

void swapNoCommDensity(double* rho, t_iterPara const * const iterPara1,
	t_iterPara const * const iterPara2, t_procData const * const procData, const int direction) {
	int fromFixedValue1, toFixedValue1, fromFixedValue2, toFixedValue2;

	if (iterPara1->fixedValue == 1) {
		fromFixedValue1 = iterPara1->fixedValue;
		toFixedValue1   = iterPara2->fixedValue+1;
		fromFixedValue2 = iterPara2->fixedValue;
		toFixedValue2   = iterPara1->fixedValue-1;
	} else {
		fromFixedValue1 = iterPara1->fixedValue;
		toFixedValue1   = iterPara2->fixedValue-1;
		fromFixedValue2 = iterPara2->fixedValue;
		toFixedValue2   = iterPara1->fixedValue+1;
	}
	for (int k = iterPara1->startOuter; k <= iterPara1->endOuter; k++) {
		for (int j = iterPara1->startInner; j <= iterPara1->endInner; j++) {

			int fromIdx1 = p_computeCellOffset(k, j, fromFixedValue1, procData->xLength, direction);
			int toIdx1   = p_computeCellOffset(k, j, toFixedValue1, procData->xLength, direction);
			int fromIdx2 = p_computeCellOffset(k, j, fromFixedValue2, procData->xLength, direction);
			int toIdx2   = p_computeCellOffset(k, j, toFixedValue2, procData->xLength, direction);

			rho[toIdx1] = rho[fromIdx1];
			rho[toIdx2] = rho[fromIdx2];
		}
	}
}


/* Function to assign iteration parameters for communication.
*
* To save communication workload:
* LEFT/RIGHT faces do only iterate the inner domain (no shared edges)
* TOP/BOTTOM faces do iterate the inner domain in X direction and the shared edges in Y direction
* FRONT/BACK faces do iterate the entire face (all shared edges)
*
* For correctness:
* If a processor lies at a real boundary the shared edges that would be normally included are excluded in this case.
*/
void p_setCommIterationParameters(t_iterPara * const iterPara, t_procData const*const procData, const int direction){

    switch(direction/2){ //integer division to get the type of face (see enum in LBDefinitions.h)
        //---------------------------------------------
        //outer = Z, inner = X, Y fixed
        //only iterate over inner domain of plane (FLUID cells)
        case 0:
        iterPara->startOuter = 1;
        iterPara->endOuter   = procData->xLength[2];
        iterPara->startInner = 1;
        iterPara->endInner   = procData->xLength[0];
        iterPara->fixedValue = (direction == LEFT) ? 1 : procData->xLength[1];
        break;

        //---------------------------------------------
        //outer = Y, inner = X, Z fixed
        //iterate over inner domain and include ghost layer in y-direction
        case 1:
        //A real boundary cell should never be included in the iteration.
        //-> Test for real boundary (MPI_PROC_NULL), and set start index accordingly
        iterPara->startOuter = 0;
        iterPara->endOuter   = procData->xLength[1]+1;
        iterPara->startInner = 1;
        iterPara->endInner   = procData->xLength[0];
        iterPara->fixedValue = (direction == BOTTOM) ? 1 : procData->xLength[2];
        break;

        //---------------------------------------------
        //outer = Z, inner = Y, X fixed
        //Iterate over entire ZY plane
        case 2:
        //A real boundary cell should never be included in the iteration.
        //-> Test for real boundary (MPI_PROC_NULL), and set start index accordingly
        iterPara->startOuter = 0;
        iterPara->endOuter   = procData->xLength[2]+1;
        iterPara->startInner = 0;
        iterPara->endInner   = procData->xLength[1]+1;
        iterPara->fixedValue = (direction == BACK) ? 1 : procData->xLength[0];
        break;

        default:
        ERROR("Invalid direction occurred. This should never happen!");
    }
}

//Function to find indices being extracted/injected
void p_assignIndices(const int face, int *const index) {

    switch (face) {
        case LEFT:     	// y = 0
        index[0] = 0;  index[1] = 5;  index[2] = 6;  index[3] = 7;  index[4] = 14; break;

        case RIGHT:     // y = xlength[1]+1
        index[0] = 4;  index[1] = 11; index[2] = 12; index[3] = 13; index[4] = 18; break;

        case TOP: 		// z = xlength[2]+1
        index[0] = 14; index[1] = 15; index[2] = 16; index[3] = 17; index[4] = 18; break;

        case BOTTOM: 	// z = 0
        index[0] = 0;  index[1] = 1;  index[2] = 2;  index[3] = 3;  index[4] = 4;  break;

        case FRONT:     // x = xlength[0]+1
        index[0] = 3;  index[1] = 7;  index[2] = 10; index[3] = 13; index[4] = 17; break;

        case BACK:      // x = 0
        index[0] = 1;  index[1] = 5;  index[2] = 8;  index[3] = 11; index[4] = 15; break;
    }
}

void finaliseMPI(t_procData const * const procData) {
    fflush(stdout);
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Rank %i: FINISHED. \n", procData->rank);
    MPI_Finalize();
}
