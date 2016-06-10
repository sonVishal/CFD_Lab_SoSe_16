#include "debug.h"
#include "helper.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>

void convertEnumWallToString(const int wall, char *wallName) {
	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	switch (wall) {
	case LEFT:
		snprintf(wallName, 80, "Proc: %d\t Wall: %s\n",myRank, "LEFT");
		break;
	case RIGHT:
		snprintf(wallName, 80, "Proc: %d\t Wall: %s\n",myRank, "RIGHT");
		break;
	case TOP:
		snprintf(wallName, 80, "Proc: %d\t Wall: %s\n",myRank, "TOP");
		break;
	case BOTTOM:
		snprintf(wallName, 80, "Proc: %d\t Wall: %s\n",myRank, "BOTTOM");
		break;
	case FRONT:
		snprintf(wallName, 80, "Proc: %d\t Wall: %s\n",myRank, "FRONT");
		break;
	case BACK:
		snprintf(wallName, 80, "Proc: %d\t Wall: %s\n",myRank, "BACK");
		break;
	default:
		snprintf(wallName, 80, "Proc: %d\t Wall: %s\n",myRank, "NULL");
		break;
	}
}


void debug_setBufferValues(double **sendBuffer, double **readBuffer, t_procData procData){
	for(int i = 0; i < 6; ++i){
		memset(sendBuffer[i], -1, procData.bufferSize[i/2]*sizeof(double));
		memset(readBuffer[i], -1, procData.bufferSize[i/2]*sizeof(double));
	}
}


#define BUFSIZE 128

int parse_output(void) {
	char command[200];
	snprintf(command, 200, "%s", "debug/ExecuteComparison.sh");

    char buf[BUFSIZE];
    FILE *fp;

    if ((fp = popen(command, "r")) == NULL) {
        printf("Error opening pipe!\n");
        return -1;
    }

    while (fgets(buf, BUFSIZE, fp) != NULL) {
        // Do whatever you want here...
        printf("OUTPUT: %s", buf);
    }

    if(pclose(fp))  {
        printf("Command not found or exited with error status\n");
        return -1;
    }

    return 1;
}


void convertEnumCellToString(const int cell, char *cellName) {
	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	switch (cell) {
	case FLUID:
		snprintf(cellName, 80, "Proc: %d\t Cell: %s\n",myRank, "FLUID");
		break;
	case NO_SLIP:
		snprintf(cellName, 80, "Proc: %d\t Cell: %s\n",myRank, "NO_SLIP");
		break;
	case MOVING_WALL:
		snprintf(cellName, 80, "Proc: %d\t Cell: %s\n",myRank, "MOVING_WALL");
		break;
	case PARALLEL_BOUNDARY:
		snprintf(cellName, 80, "Proc: %d\t Cell: %s\n",myRank, "PARALLEL_BOUNDARY");
		break;
	default:
		snprintf(cellName, 80, "Proc: %d\t Cell: %s\n",myRank, "NULL");
		break;
	}
}

void printWallEnum(const int wall) {
	char wallName[80];
	convertEnumWallToString(wall,wallName);
	fprintf(stderr, "%s\n",wallName);
	fflush(stdout);
	fflush(stderr);
}

void printCellEnum(const int cellType) {
	char typeName[80];
	convertEnumWallToString(cellType,typeName);
	fprintf(stderr, "%s\n",typeName);
	fflush(stdout);
	fflush(stderr);
}

void printProcData(t_procData procData) {
	fprintf(stderr,"------------- Proc %d -------------\n",procData.rank);
	fprintf(stderr,"My length = (%d,%d,%d)\n",procData.xLength[0],procData.xLength[1],procData.xLength[2]);
	fprintf(stderr,"My neighbors are\n");
	fprintf(stderr,"LEFT = %d, RIGHT = %d,\nTOP = %d, BOTTOM = %d,\nFRONT = %d, BACK = %d\n",procData.neighbours[LEFT],procData.neighbours[RIGHT],procData.neighbours[TOP],procData.neighbours[BOTTOM],procData.neighbours[FRONT],procData.neighbours[BACK]);
	fprintf(stderr,"-----------------------------------\n");
	fflush(stdout);
	fflush(stderr);
}

void printProcDataPos(t_procData procData, int *pos) {
	fprintf(stderr,"------------- Proc %d -------------\n",procData.rank);
	fprintf(stderr,"My position = (%d,%d,%d)\n",pos[0],pos[1],pos[2]);
	fprintf(stderr,"My length = (%d,%d,%d)\n",procData.xLength[0],procData.xLength[1],procData.xLength[2]);
	fprintf(stderr,"My neighbors are\n");
	fprintf(stderr,"LEFT\t= %d,\tRIGHT\t= %d,\nTOP\t= %d,\tBOTTOM\t= %d,\nFRONT\t= %d,\tBACK\t= %d\n",procData.neighbours[LEFT],procData.neighbours[RIGHT],procData.neighbours[TOP],procData.neighbours[BOTTOM],procData.neighbours[FRONT],procData.neighbours[BACK]);
	fprintf(stderr,"-----------------------------------\n");
	fflush(stdout);
	fflush(stderr);
}

void writeVtkOutputDebug(const double * const collideField,
		const int * const flagField, const char * filename,
		unsigned int t, t_procData procData, int *procsPerAxis)
{
	// Files related variables
	char pFileName[80];
	FILE *fp = NULL;

	// Create the file with time information in the name
	sprintf(pFileName, "%s.%i.%i.vts", "debug/Debug",procData.rank,t);
	fp  = fopen(pFileName, "w");

	int myPos[3] = {0,0,0};
	// printf("Rank %d\n",procData.rank);
	p_rankToPos(procsPerAxis, procData.rank, myPos);

	// Check if files were opened or not
	if(fp == NULL)
	{
		char szBuff[80];
		sprintf(szBuff, "Failed to open %s", pFileName);
		ERROR(szBuff);
		return;
	}

	fprintf(fp,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");

	// Write the point data for the domain
	writevtkPointCoordinatesDebug(fp,procData.xLength,myPos);

	fprintf(fp,"<CellData>\n");
	fprintf(fp,"<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"BoundaryType\">\n");

	int x, y, z;            // iteration variables
	int idx;                // cell index

	// Temporary variables for (xlength+2)^2
	int const xylen = (procData.xLength[0]+2)*(procData.xLength[1]+2);

	// Temporary variables for z and y offsets
	int zOffset, yzOffset;
	// Write cell average denisties the file
	for(z = 0; z <= procData.xLength[2]+1; z++) {
		zOffset = z*xylen;
		for(y = 0; y <= procData.xLength[1]+1; y++) {
			yzOffset = zOffset + y*(procData.xLength[0]+2);
			for(x = 0; x <= procData.xLength[0]+1; x++) {
				// Compute the base index for collideField
				idx = (yzOffset + x);

				// Write Boundary type
				fprintf(fp, "%d\n", flagField[idx]);
			}
		}
	}
	fprintf(fp,"</DataArray>\n");
	fprintf(fp,"</CellData>\n");
	fprintf(fp,"</Piece>\n");
	fprintf(fp,"</StructuredGrid>\n");
	fprintf(fp,"</VTKFile>\n");
	// Close file
	if(fclose(fp))
	{
		char szBuff[80];
		sprintf(szBuff, "Failed to close %s", pFileName);
		ERROR(szBuff);
	}
}

void writevtkPointCoordinatesDebug(FILE *fp, int *xlength, int *myPos) {
	int x, y, z;
	// printf("Position = (%d,%d,%d)\n",myPos[0],myPos[1],myPos[2]);
	unsigned int x1 = myPos[0]*(xlength[0]+2);
	unsigned int x2 = x1 + xlength[0] + 2;
	unsigned int y1 = myPos[1]*(xlength[1]+2);
	unsigned int y2 = y1 + xlength[1] + 2;
	unsigned int z1 = myPos[2]*(xlength[2]+2);
	unsigned int z2 = z1 + xlength[2] + 2;
	fprintf(fp,"<StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n",x1,x2,y1,y2,z1,z2);
	fprintf(fp,"<Piece Extent=\"%d %d %d %d %d %d\">\n",x1,x2,y1,y2,z1,z2);
	fprintf(fp,"<Points>\n");
	fprintf(fp,"<DataArray type=\"UInt32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	// We have xlength + 1 points for xlength cells in each direction
	for(z = z1; z <= z2; z++) {
		for(y = y1; y <= y2; y++) {
			for(x = x1; x <= x2; x++) {
				fprintf(fp, "%d %d %d\n", x, y, z);
			}
		}
	}
	fprintf(fp,"</DataArray>\n");
	fprintf(fp,"</Points>\n");
}

void p_writeCombinedPVTSFileDebug(const char * filename, unsigned int t, int xlength, int *procsPerAxis) {
	// Files related variables
	char pFileName[80];
	FILE *fp = NULL;

	// Create the file with time information in the name
	sprintf(pFileName, "%s.%i.pvts", "debug/Debug",t);
	fp  = fopen(pFileName, "w");

	// Check if files were opened or not
	if(fp == NULL)
	{
		char szBuff[80];
		sprintf(szBuff, "Failed to open %s", pFileName);
		ERROR(szBuff);
		return;
	}

	fprintf(fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(fp, "<PStructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\" GhostLevel=\"1\">\n",xlength+4,xlength+4,xlength+4);
	fprintf(fp, "<PPoints>\n");
	fprintf(fp, "%s\n","<PDataArray NumberOfComponents=\"3\" type=\"UInt32\" />");
	fprintf(fp, "</PPoints>\n");
	fprintf(fp, "<PCellData>\n");
	fprintf(fp, "<PDataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"BoundaryType\"/>\n");
	fprintf(fp, "</PCellData>\n");

	// Perform domain decomposition again
	int procXlength[3] = {0,0,0};
	procXlength[0] = (xlength+4)/procsPerAxis[0];
	procXlength[1] = (xlength+4)/procsPerAxis[1];
	procXlength[2] = (xlength+4)/procsPerAxis[2];
	int x1,x2,y1,y2,z1,z2;

	for (int k = 0; k < procsPerAxis[2]; k++) {
		for (int j = 0; j < procsPerAxis[1]; j++) {
			for (int i = 0; i < procsPerAxis[0]; i++) {
				procXlength[0] += ((procsPerAxis[0]-1)==i)?(xlength+4)%procsPerAxis[0]:0;
				procXlength[1] += ((procsPerAxis[1]-1)==j)?(xlength+4)%procsPerAxis[1]:0;
				procXlength[2] += ((procsPerAxis[2]-1)==k)?(xlength+4)%procsPerAxis[2]:0;
				x1 = i*procXlength[0]; x2 = x1 + procXlength[0];
				y1 = j*procXlength[1]; y2 = y1 + procXlength[1];
				z1 = k*procXlength[2]; z2 = z1 + procXlength[2];
				snprintf(pFileName, 80, "%s.%i.%i.vts",filename,i+(j+k*procsPerAxis[1])*procsPerAxis[0],t);
				// This is stupid "../<fileName>" but what the heck
				fprintf(fp, "<Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s.%i.%i.vts\"/>\n",x1,x2,y1,y2,z1,z2,"Debug",i+j*procsPerAxis[0]+k*procsPerAxis[0]*procsPerAxis[1],t);
			}
		}
	}
	fprintf(fp, "</PStructuredGrid>\n");
	fprintf(fp, "</VTKFile>\n");

	// Close file
	if(fclose(fp))
	{
		char szBuff[80];
		sprintf(szBuff, "Failed to close %s", pFileName);
		ERROR(szBuff);
	}
}


void checkCollideFieldDebug(char *referenceFile, double *currentCollideField, int size){
	FILE *fp = NULL;

	char filenameExt[100];
	snprintf(filenameExt, 100, "%s.array", referenceFile);

	fp  = fopen(filenameExt, "r");
	if(fp == NULL)
	{
		char szBuff[200];
		sprintf(szBuff, "Failed to open %s. \nMaybe you have to run "
				"'writeCollideFieldDebug' with the debug reference scenario first.", filenameExt);
		ERROR(szBuff);
		return;
	}

	int counter = 0;
	const double TOL = 1e-10;

	while(1){
		double refVal, curVal; //have to correspond!
		int ret = fscanf(fp, "%lf", &refVal);

		if(ret == EOF || counter > size){
			//        	printf("ret=%i, counter=%i, size=%i", ret, counter, size);
			if(counter == size-1 && ret == EOF) 	break; //no error
			else  								ERROR("TODO");
		}

		curVal = currentCollideField[counter];
		//        printf("TOL=%f, counter=%i refVal: %.16f, curVal: %.16f", TOL, counter, refVal, curVal);
		// Check for relative error to get number of correct digits
		// rather than absolute error
		if(fabs((refVal - curVal)/refVal) > TOL){
			char msg[200];
			snprintf(msg, 200, "TOL=%f, counter=%i refVal: %.16f, curVal: %.16f", TOL, counter, refVal, curVal);
			ERROR(msg);
		}
		counter++;
	}

	if(fclose(fp))
	{
		char szBuff[80];
		sprintf(szBuff, "Failed to close %s", filenameExt);
		ERROR(szBuff);
	}

	printf("INFO: CHECK SUCCESSFUL \n");
}
