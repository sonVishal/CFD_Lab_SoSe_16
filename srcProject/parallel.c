#include "parallel.h"

void initialiseMPI(int *rank, int *numRanks, int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,numRanks);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);
}

/*
* Divides the entire domain into sub-domains. The function sets the processor specific data to determine the size
* of the sub-domain and determines the neighbours.
*/
void domainDecompositionAndNeighbors(t_procData *const procData, const int xlength, int const * const procsPerAxis) {

	#define M_NBGH procData->neighbours
	#define M_PDICNBGH procData->periodicNeighbours
	#define M_EDGENBGH procData->periodicEdgeNeighbours
	#define M_ISEDGE(N1, N2) M_NBGH[N1] == MPI_PROC_NULL && M_NBGH[N2] == MPI_PROC_NULL

	int procPos[3] = {0,0,0};
	// Get the position of the process based on its rank
	p_rankToPos(procsPerAxis,procData->rank,procPos);

	/* Compute the subdomain size and save it into procData */
	int baseLength[3] = {xlength/procsPerAxis[0], xlength/procsPerAxis[1], xlength/procsPerAxis[2]};

	// If the proc is at the end of some axis then add the remaining length
	procData->xLength[0] = (procPos[0] == procsPerAxis[0]-1)?xlength-(procsPerAxis[0]-1)*baseLength[0]:baseLength[0];
	procData->xLength[1] = (procPos[1] == procsPerAxis[1]-1)?xlength-(procsPerAxis[1]-1)*baseLength[1]:baseLength[1];
	procData->xLength[2] = (procPos[2] == procsPerAxis[2]-1)?xlength-(procsPerAxis[2]-1)*baseLength[2]:baseLength[2];

	/* Set neighbours for parallel boundaries */
	/* Decide whether it is a ghost boundary (= MPI_PROC_NULL) or a parallel boundary (assign rank of neighbor) */
	M_NBGH[LEFT]   = (procPos[1] == 0) 				   ? MPI_PROC_NULL : procData->rank-procsPerAxis[0];
	M_NBGH[RIGHT]  = (procPos[1] == procsPerAxis[1]-1) ? MPI_PROC_NULL : procData->rank+procsPerAxis[0];

	M_NBGH[BACK]   = (procPos[0] == 0) 				   ? MPI_PROC_NULL : procData->rank-1;
	M_NBGH[FRONT]  = (procPos[0] == procsPerAxis[0]-1) ? MPI_PROC_NULL : procData->rank+1;

	M_NBGH[BOTTOM] = (procPos[2] == 0) 				   ? MPI_PROC_NULL : procData->rank-procsPerAxis[1]*procsPerAxis[0];
	M_NBGH[TOP]    = (procPos[2] == procsPerAxis[2]-1) ? MPI_PROC_NULL : procData->rank+procsPerAxis[1]*procsPerAxis[0];

	/* Set neighbours for periodic boundaries */
	M_PDICNBGH[LEFT]   = M_NBGH[LEFT]==MPI_PROC_NULL  ? procData->rank+procsPerAxis[0]*(procsPerAxis[1]-1) : MPI_PROC_NULL;
	M_PDICNBGH[RIGHT]  = M_NBGH[RIGHT]==MPI_PROC_NULL ? procData->rank-procsPerAxis[0]*(procsPerAxis[1]-1) : MPI_PROC_NULL;

	M_PDICNBGH[BACK]   = M_NBGH[BACK]==MPI_PROC_NULL  ? procData->rank-procsPerAxis[0]+1 : MPI_PROC_NULL;
	M_PDICNBGH[FRONT]  = M_NBGH[FRONT]==MPI_PROC_NULL ? procData->rank+procsPerAxis[0]-1 : MPI_PROC_NULL;

	M_PDICNBGH[BOTTOM] = M_NBGH[BOTTOM]==MPI_PROC_NULL ? procData->rank+procsPerAxis[0]*procsPerAxis[1]*(procsPerAxis[2]-1) : MPI_PROC_NULL;
	M_PDICNBGH[TOP]    = M_NBGH[TOP]==MPI_PROC_NULL    ? procData->rank-procsPerAxis[0]*procsPerAxis[1]*(procsPerAxis[2]-1) : MPI_PROC_NULL;

	/* Set neighbours for shared edges */
	M_EDGENBGH[0] = M_ISEDGE(LEFT, BOTTOM)  ? M_PDICNBGH[TOP]+procsPerAxis[0]*(procsPerAxis[1]-1): MPI_PROC_NULL;
	M_EDGENBGH[1] = M_ISEDGE(FRONT, BOTTOM) ? M_PDICNBGH[TOP]-(procsPerAxis[0]-1): MPI_PROC_NULL;
	M_EDGENBGH[2] = M_ISEDGE(RIGHT, BOTTOM) ? M_PDICNBGH[TOP]-(procsPerAxis[0]*(procsPerAxis[1]-1)): MPI_PROC_NULL;
	M_EDGENBGH[3] = M_ISEDGE(BACK, BOTTOM)  ? M_PDICNBGH[TOP]+(procsPerAxis[0]-1): MPI_PROC_NULL;

	M_EDGENBGH[4] = M_ISEDGE(BACK, LEFT)    ? M_PDICNBGH[LEFT]+(procsPerAxis[0]-1): MPI_PROC_NULL;
	M_EDGENBGH[5] = M_ISEDGE(LEFT, FRONT)   ? M_PDICNBGH[LEFT]-(procsPerAxis[0]-1): MPI_PROC_NULL;
	M_EDGENBGH[6] = M_ISEDGE(FRONT, RIGHT)  ? M_PDICNBGH[RIGHT]-(procsPerAxis[0]-1): MPI_PROC_NULL;
	M_EDGENBGH[7] = M_ISEDGE(RIGHT, BACK)   ? M_PDICNBGH[RIGHT]+(procsPerAxis[0]-1): MPI_PROC_NULL;

	M_EDGENBGH[8] = M_ISEDGE(TOP, LEFT)     ? M_PDICNBGH[BOTTOM]+(procsPerAxis[0]*(procsPerAxis[1]-1)): MPI_PROC_NULL;
	M_EDGENBGH[9] = M_ISEDGE(TOP, FRONT)    ? M_PDICNBGH[BOTTOM]-(procsPerAxis[0]-1): MPI_PROC_NULL;
	M_EDGENBGH[10] = M_ISEDGE(TOP, RIGHT)   ? M_PDICNBGH[BOTTOM]-(procsPerAxis[0]*(procsPerAxis[1]-1)): MPI_PROC_NULL;
	M_EDGENBGH[11] = M_ISEDGE(TOP, BACK)    ? M_PDICNBGH[BOTTOM]+(procsPerAxis[0]-1): MPI_PROC_NULL;
}

/*
* Initialize sendBuffer and readBuffer. They are only initialized in a certain direction when there is a valid neighbor
* i.e. no domain boundary
*/
void initialiseBuffers(double *sendBuffer[6], double *readBuffer[6], int const * const xlength,
	int const * const neighbours, int procBufferSize[3]) {

    const int nrDistSwap = 5;
    const int db = sizeof(double); //double in bytes

    /*TODO: (DL) Do we even need for each side a new buffer?
    When we only have one readBuffer and one sendBuffer, and we allocate it such that also
    the shared edges are included, then we can reuse this buffer for all directions?
    */
    
    // XZ inner domain (no edges included)
    int bufferSize		= nrDistSwap*(xlength[0]*xlength[2]);
    procBufferSize[0] 	= bufferSize; //Valid for left and right
    sendBuffer[LEFT] 	= (neighbours[LEFT]  != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;
    sendBuffer[RIGHT] 	= (neighbours[RIGHT] != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;
    readBuffer[LEFT] 	= (neighbours[LEFT]  != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;
    readBuffer[RIGHT] 	= (neighbours[RIGHT] != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;

    // XY plane including edges at the boundary to left/right
    bufferSize 			= nrDistSwap*(xlength[0]*(xlength[1]+2));
    procBufferSize[1] 	= bufferSize; //Valid for top and bottom
    sendBuffer[TOP] 	= (neighbours[TOP] 	  != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;
    sendBuffer[BOTTOM] 	= (neighbours[BOTTOM] != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;
    readBuffer[TOP] 	= (neighbours[TOP]    != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;
    readBuffer[BOTTOM] 	= (neighbours[BOTTOM] != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;

    // YZ plane including all edges
    bufferSize 			= nrDistSwap*((xlength[1]+2)*(xlength[2]+2));
    procBufferSize[2] 	= bufferSize; //Valid for front and back
    sendBuffer[FRONT] 	= (neighbours[FRONT]  != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;
    sendBuffer[BACK] 	= (neighbours[BACK]   != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;
    readBuffer[FRONT] 	= (neighbours[FRONT]  != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;
    readBuffer[BACK] 	= (neighbours[BACK]   != MPI_PROC_NULL)?(double *) malloc(bufferSize*db):NULL;
}

/*
* Main algorithm to carry out the communication between the processes.
*/
void communicate(double** sendBuffer, double**readBuffer, double* collideField, t_procData const * const procData){
    //Run extract, swap, inject for all sides and cells.

    int 		index1[5], index2[5];
    t_iterPara  iterPara1, iterPara2;

    // Iterate through directions (left/right, top/bottom, front/back)
    for (int direction = LEFT; direction <= BACK; direction+=2) {

        const int face1 = direction;   //LEFT,  TOP,    FRONT
        const int face2 = direction+1; //RIGHT, BOTTOM, BACK

        /* Set iteration parameters for both faces of a direction */
        p_setCommIterationParameters(&iterPara1, procData, face1);
        p_setCommIterationParameters(&iterPara2, procData, face2);

        /* Assign the indices for the 5 distributions that go out of the two faces of the current direction */
        p_assignIndices(face1, index1);
        p_assignIndices(face2, index2);

        /* Extract the distributions from collideField to the send buffer */
        if(procData->neighbours[face1] != MPI_PROC_NULL)
            extract(sendBuffer[face1], collideField, &iterPara1, procData, face1, index1);

        if(procData->neighbours[face2] != MPI_PROC_NULL)
            extract(sendBuffer[face2], collideField, &iterPara2, procData, face2, index2);

        /* Swap the distributions */
        swap(sendBuffer, readBuffer, procData, direction);

        /* Inject the distributions from read buffer to collide field */
        //NOTE: opposite index get injected from the readBuffer
        if(procData->neighbours[face1] != MPI_PROC_NULL)
            inject(readBuffer[face1], collideField, &iterPara1, procData, face1, index2);

        if(procData->neighbours[face2] != MPI_PROC_NULL)
            inject(readBuffer[face2], collideField, &iterPara2, procData, face2, index1);
    }
}

/*
* Copy distributions into sendBuffer.
*/
void extract(double sendBuffer[], double const * const collideField, t_iterPara const * const iterPara, t_procData const * const procData,
    const int direction, int const * const index){

    int currentIndexField = -1;//initially invalid assignment
    int currentIndexBuff  = 0; //simple counter

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
                sendBuffer[currentIndexBuff++] = collideField[currentIndexField+index[i]];
            }
        }
    }
}

/*
* Swap distributions with neighboring processes.
*/
void swap(double** sendBuffer, double** readBuffer, t_procData const * const procData, const int direction){

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
    if(procData->neighbours[direction] != MPI_PROC_NULL){
        MPI_Sendrecv(sendBuffer[direction], bufferSize , MPI_DOUBLE, proc1, 0, readBuffer[direction],
            bufferSize, MPI_DOUBLE, proc1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //Send proc2 receive proc2
    if(procData->neighbours[direction+1] != MPI_PROC_NULL){
        MPI_Sendrecv(sendBuffer[direction+1], bufferSize , MPI_DOUBLE, proc2, 0, readBuffer[direction+1],
            bufferSize, MPI_DOUBLE, proc2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

/*
* Copy read buffer into ghost layer of collideField.
*/
void inject(double const * const readBuffer, double* collideField, t_iterPara *const iterPara, t_procData const * const procData,
    const int direction, int const * const index){

    int currentIndexField = -1; //initially assign to invalid
    int currentIndexBuff = 0;

    //    switch (direction){
    //        case LEFT:   iterPara->fixedValue--; break;
    //        case RIGHT:  iterPara->fixedValue++; break;
    //        case TOP:    iterPara->fixedValue++; break;
    //        case BOTTOM: iterPara->fixedValue--; break;
    //        case FRONT:  iterPara->fixedValue++; break;
    //        case BACK:   iterPara->fixedValue--; break;
    //    }
    static int shiftFixedValue[6] = {-1, 1, 1, -1, 1, -1}; //static to allocate/assign only once
    iterPara->fixedValue += shiftFixedValue[direction]; //change value also in outer scope okay, because not needed afterwards

    //k - corresponds to the 'outer' value when computing the offset
    for(int k = iterPara->startOuter; k <= iterPara->endOuter; ++k){
        //j - corresponds to the 'inner' value
        for(int j = iterPara->startInner; j <= iterPara->endInner; ++j){

            currentIndexField  = Q*p_computeCellOffset(k, j, iterPara->fixedValue, procData->xLength, direction);

            //out of bound checks
            assert(currentIndexField < Q*(procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2)
            && currentIndexField >= 0);

            for (int i = 0; i < 5; i++) {
                //out of bound checks
                assert(currentIndexBuff < procData->bufferSize[direction/2] && currentIndexBuff >=0);
                collideField[currentIndexField+index[i]] = readBuffer[currentIndexBuff++];
            }
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
        iterPara->startOuter = (procData->neighbours[LEFT]  == MPI_PROC_NULL) ? 1:0;
        iterPara->endOuter   = (procData->neighbours[RIGHT] == MPI_PROC_NULL) ? procData->xLength[1] : procData->xLength[1]+1;
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
        iterPara->startOuter = (procData->neighbours[BOTTOM] == MPI_PROC_NULL) ? 1:0;
        iterPara->endOuter   = (procData->neighbours[TOP]    == MPI_PROC_NULL) ? procData->xLength[2] : procData->xLength[2]+1;
        iterPara->startInner = (procData->neighbours[LEFT]   == MPI_PROC_NULL) ? 1:0;
        iterPara->endInner   = (procData->neighbours[RIGHT]  == MPI_PROC_NULL) ? procData->xLength[1] : procData->xLength[1]+1;
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
        index[0] = 0; index[1] = 5; index[2] = 6; index[3] = 7; index[4] = 14;
        break;

        case RIGHT:     // y = xlength[1]+1
        index[0] = 4; index[1] = 11; index[2] = 12; index[3] = 13; index[4] = 18;
        break;

        case TOP: 		// z = xlength[2]+1
        index[0] = 14; index[1] = 15; index[2] = 16; index[3] = 17; index[4] = 18;
        break;

        case BOTTOM: 	// z = 0
        index[0] = 0; index[1] = 1; index[2] = 2; index[3] = 3; index[4] = 4;
        break;

        case FRONT:     // x = xlength[0]+1
        index[0] = 3; index[1] = 7; index[2] = 10; index[3] = 13; index[4] = 17;
        break;

        case BACK:      // x = 0
        index[0] = 1; index[1] = 5; index[2] = 8; index[3] =  11; index[4] = 15;
        break;
    }
}

void finaliseMPI(t_procData const * const procData) {
    fflush(stdout);
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Rank %i: FINISHED. \n", procData->rank);
    MPI_Finalize();
}
