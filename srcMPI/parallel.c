#include "parallel.h"

void initialiseMPI(int *rank, int *numRanks, int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,numRanks);
    MPI_Comm_rank(MPI_COMM_WORLD,rank);
}

void p_domainDecompositionAndNeighbors(t_procData *procData, const int xlength, const int * const procsPerAxis) {

	int procPos[3] = {0,0,0};
	// Get the position of the process based on its rank
	p_rankToPos(procsPerAxis,procData->rank,procPos);

    /* Compute the subdomain size and save it into procData */
 	procData->xLength[0] = xlength/procsPerAxis[0];
    procData->xLength[1] = xlength/procsPerAxis[1];
    procData->xLength[2] = xlength/procsPerAxis[2];
    // If the proc is at the end of some axis then add the remaining length
    procData->xLength[0] += (procPos[0] == procsPerAxis[0]-1)?xlength%procsPerAxis[0]:0;
    procData->xLength[1] += (procPos[1] == procsPerAxis[1]-1)?xlength%procsPerAxis[1]:0;
    procData->xLength[2] += (procPos[2] == procsPerAxis[2]-1)?xlength%procsPerAxis[2]:0;

    /* Decide whether it is a ghost boundary (= MPI_PROC_NULL) or a parallel boundary (rank of neighbour) */
    procData->neighbours[LEFT]   = (procPos[1] == 0) 				  ? MPI_PROC_NULL : procData->rank-procsPerAxis[0];
    procData->neighbours[RIGHT]  = (procPos[1] == procsPerAxis[1]-1)  ? MPI_PROC_NULL : procData->rank+procsPerAxis[0];

    procData->neighbours[BACK]   = (procPos[0] == 0) 				  ? MPI_PROC_NULL : procData->rank-1;
    procData->neighbours[FRONT]  = (procPos[0] == procsPerAxis[0]-1)  ? MPI_PROC_NULL : procData->rank+1;

    procData->neighbours[BOTTOM] = (procPos[2] == 0) 				  ? MPI_PROC_NULL : procData->rank-procsPerAxis[1]*procsPerAxis[0];
    procData->neighbours[TOP]    = (procPos[2] == procsPerAxis[2]-1)  ? MPI_PROC_NULL : procData->rank+procsPerAxis[1]*procsPerAxis[0];
}

void initialiseBuffers(double *sendBuffer[6], double *readBuffer[6], int *xlength, int *neighbours, int procBufferSize[3]) {

	const int nrDistSwap = 5;
	const int db = sizeof(double); //double in bytes

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

void communicate(double** sendBuffer, double**readBuffer, double* collideField, const t_procData *procData){
    //Run extract, swap, inject for all sides and cells.
    //
    int index1[5];
    int index2[5];

    t_iterPara iterPara1;
    t_iterPara iterPara2;

	// Treat opposite directions in one go
    for (int direction = LEFT; direction <= BACK; direction+=2) {

		// Set iteration parameters for the opposite directions
        p_setCommIterationParameters(&iterPara1, procData, direction);
        p_setCommIterationParameters(&iterPara2, procData, direction+1);

		// Assign the indices for the 5 distributions that go out of the cell in the direction
        p_assignIndices(direction,   index1);
        p_assignIndices(direction+1, index2);

		// Extract the distributions from collide field to the send buffer
        if(procData->neighbours[direction] != MPI_PROC_NULL)
            extract(sendBuffer[direction], collideField, &iterPara1, procData, direction, index1);

        if(procData->neighbours[direction+1] != MPI_PROC_NULL)
            extract(sendBuffer[direction+1], collideField, &iterPara2, procData, direction+1, index2);

		// Swap the distributions
        swap(sendBuffer, readBuffer, procData, direction);

		// Inject the distributions from read buffer to collide field
        if(procData->neighbours[direction] != MPI_PROC_NULL)
            inject(readBuffer[direction], collideField, &iterPara1, procData, direction, index2);

        if(procData->neighbours[direction+1] != MPI_PROC_NULL)
            inject(readBuffer[direction+1], collideField, &iterPara2, procData, direction+1, index1);
    }

}

//Copy distributions needed to sendbuffer.
void extract( double sendBuffer[], double* collideField, const t_iterPara * const iterPara, const t_procData *procData,
              int direction, int* index){

    int currentIndexField;
    int currentIndexBuff= 0;

#ifndef NDEBUG
    int fieldSize = Q*(procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2);
#endif

    //k - corresponds to the 'outer' value when computing the offset
    for(int k = iterPara->startOuter; k <= iterPara->endOuter; ++k){
        //j - corresponds to the 'inner' value
        for(int j = iterPara->startInner; j <= iterPara->endInner; ++j){

            currentIndexField  = Q*p_computeCellOffset(k, j, iterPara->fixedValue, procData->xLength, direction);

            assert(currentIndexBuff < procData->bufferSize[direction/2] && currentIndexBuff >=0);
            assert(currentIndexField < fieldSize  && currentIndexField >= 0);
            for (int i = 0; i < 5; i++) {
                sendBuffer[currentIndexBuff++] = collideField[currentIndexField+index[i]];
            }
        }
    }

}

//Send distributions and wait to receive.
void swap(double** sendBuffer, double** readBuffer, const t_procData *procData, int direction){
//int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //int dest, int sendtag, void *recvbuf, int recvcount,
    //MPI_Datatype recvtype, int source, int recvtag,
    //MPI_Comm comm, MPI_Status *status)

    //temp variables for readabillity.
    int bufferSize = procData->bufferSize[direction/2];
    int proc1 = procData->neighbours[direction];
    int proc2 = procData->neighbours[direction+1];

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

//Copy read buffer to ghost layer
void inject(double readBuffer[], double* collideField, t_iterPara *iterPara, const t_procData *procData,
            int direction, int *index){

    int currentIndexField;
    int currentIndexBuff= 0;

#ifndef NDEBUG
    int fieldSize = Q*(procData->xLength[0]+2)*(procData->xLength[1]+2)*(procData->xLength[2]+2);
#endif

    int shiftFixedValue[6] = {-1, 1, 1, -1, 1, -1};
    iterPara->fixedValue += shiftFixedValue[direction];

//    switch (direction){
//        case LEFT:
//              iterPara->fixedValue--;
//              break;
//        case RIGHT:
//              iterPara->fixedValue++;
//              break;
//        case TOP:
//              iterPara->fixedValue++;
//              break;
//        case BOTTOM:
//              iterPara->fixedValue--;
//              break;
//        case FRONT:
//              iterPara->fixedValue++;
//              break;
//        case BACK:
//              iterPara->fixedValue--;
//              break;
//    }

    //k - corresponds to the 'outer' value when computing the offset
    for(int k = iterPara->startOuter; k <= iterPara->endOuter; ++k){
        //j - corresponds to the 'inner' value
        for(int j = iterPara->startInner; j <= iterPara->endInner; ++j){

            currentIndexField  = Q*p_computeCellOffset(k, j, iterPara->fixedValue, procData->xLength, direction);

            assert(currentIndexBuff < procData->bufferSize[direction/2] && currentIndexBuff >=0);
            assert(currentIndexField < fieldSize  && currentIndexField >= 0);
            for (int i = 0; i < 5; i++) {
                collideField[currentIndexField+index[i]] = readBuffer[currentIndexBuff++];
            }
        }
    }
}

//Function to assign iteration parameters for communication.
void p_setCommIterationParameters(t_iterPara *iterPara, const t_procData *procData, const int direction){

	switch(direction){
	//---------------------------------------------
	//outer = Z, inner = X, Y fixed
    //only iterate over inner domain of plane (FLUID cells)
	case LEFT:
        iterPara->startOuter = 1;
		iterPara->endOuter   = procData->xLength[2];
        iterPara->startInner = 1;
		iterPara->endInner   = procData->xLength[0];
		iterPara->fixedValue = 1;
		break;
	case RIGHT:
        iterPara->startOuter = 1;
		iterPara->endOuter   = procData->xLength[2];
        iterPara->startInner = 1;
		iterPara->endInner   = procData->xLength[0];
		iterPara->fixedValue = procData->xLength[1];
		break;
	//---------------------------------------------
	//outer = Y, inner = X, Z fixed
    //iterate over inner domain and include ghost layer in y-direction
	case TOP:
        iterPara->startOuter = 0;
		iterPara->endOuter   = procData->xLength[1]+1;
        iterPara->startInner = 1;
		iterPara->endInner   = procData->xLength[0];
		iterPara->fixedValue = procData->xLength[2];
		break;
	case BOTTOM:
        iterPara->startOuter = 0;
		iterPara->endOuter   = procData->xLength[1]+1;
        iterPara->startInner = 1;
		iterPara->endInner   = procData->xLength[0];
		iterPara->fixedValue = 1;
		break;

	//---------------------------------------------
	//outer = Z, inner = Y, X fixed
    //Iterate over entire ZY plane
	case FRONT:
        iterPara->startOuter = 0;
		iterPara->endOuter   = procData->xLength[2]+1;
        iterPara->startInner = 0;
		iterPara->endInner   = procData->xLength[1]+1;
		iterPara->fixedValue = procData->xLength[0];
		break;
	case BACK:
        iterPara->startOuter = 0;
		iterPara->endOuter = procData->xLength[2]+1;
        iterPara->startInner = 0;
		iterPara->endInner = procData->xLength[1]+1;
		iterPara->fixedValue = 1;
		break;

	default:
		ERROR("Invalid direction occurred. This should never happen!");
	}
}

//Function to find indecies being extracted/injected
void p_assignIndices(int direction, int *index) {
	switch (direction) {
		case TOP:
			// z = xlength[2]+1
			index[0] = 14; index[1] = 15; index[2] = 16; index[3] = 17; index[4] = 18;
			break;
		case BOTTOM:
			// z = 0
			index[0] = 0; index[1] = 1; index[2] = 2; index[3] = 3; index[4] = 4;
			break;
		case FRONT:
			// x = xlength[0]+1
			index[0] = 3; index[1] = 7; index[2] = 10; index[3] = 13; index[4] = 17;
			break;
		case BACK:
			// x = 0
			index[0] = 1; index[1] = 5; index[2] = 8; index[3] =  11; index[4] = 15;
			break;
		case LEFT:
			// y = 0
			index[0] = 0; index[1] = 5; index[2] = 6; index[3] = 7; index[4] = 14;
			break;
		case RIGHT:
			// y = xlength[1]+1
			index[0] = 4; index[1] = 11; index[2] = 12; index[3] = 13; index[4] = 18;
			break;
	}
}

void finaliseMPI(t_procData *procData) {
    fflush(stdout);
    fflush(stderr);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Rank %i: FINISHED. \n", procData->rank);
    MPI_Finalize();
}

void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}
