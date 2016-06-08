#include "debug.h"
#include "helper.h"
#include <stdio.h>

void printProcData(t_procData procData) {
    printf("------------- Proc %d -------------\n",procData.rank);
    printf("My length = (%d,%d,%d)\n",procData.xLength[0],procData.xLength[1],procData.xLength[2]);
    printf("My neighbors are\n");
    printf("LEFT = %d, RIGHT = %d,\nTOP = %d, BOTTOM = %d,\nFRONT = %d, BACK = %d\n",procData.neighbours[LEFT],procData.neighbours[RIGHT],procData.neighbours[TOP],procData.neighbours[BOTTOM],procData.neighbours[FRONT],procData.neighbours[BACK]);
    printf("-----------------------------------\n");
}

void printProcDataPos(t_procData procData, int *pos) {
    printf("------------- Proc %d -------------\n",procData.rank);
    printf("My position = (%d,%d,%d)\n",pos[0],pos[1],pos[2]);
    printf("My length = (%d,%d,%d)\n",procData.xLength[0],procData.xLength[1],procData.xLength[2]);
    printf("My neighbors are\n");
    printf("LEFT\t= %d,\tRIGHT\t= %d,\nTOP\t= %d,\tBOTTOM\t= %d,\nFRONT\t= %d,\tBACK\t= %d\n",procData.neighbours[LEFT],procData.neighbours[RIGHT],procData.neighbours[TOP],procData.neighbours[BOTTOM],procData.neighbours[FRONT],procData.neighbours[BACK]);
    printf("-----------------------------------\n");
}

void writeVtkOutputDebug(const double * const collideField,
    const int * const flagField, const char * filename,
    unsigned int t, t_procData procData, int *procsPerAxis)
{
    // Files related variables
    char pFileName[80];
    FILE *fp = NULL;

    // Create the file with time information in the name
    sprintf(pFileName, "%s.%i.%i.vts", "pv_files/Debug",procData.rank,t);
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
    sprintf(pFileName, "%s.%i.pvts", "pv_files/Debug",t);
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
    procXlength[0] = xlength/procsPerAxis[0];
    procXlength[1] = xlength/procsPerAxis[1];
    procXlength[2] = xlength/procsPerAxis[2];
    int x1,x2,y1,y2,z1,z2;

    for (int k = 0; k < procsPerAxis[2]; k++) {
        for (int j = 0; j < procsPerAxis[1]; j++) {
            for (int i = 0; i < procsPerAxis[0]; i++) {
                procXlength[0] += ((procsPerAxis[0]-1)==i)?(xlength+2)%procsPerAxis[0]:0;
                procXlength[1] += ((procsPerAxis[1]-1)==j)?(xlength+2)%procsPerAxis[1]:0;
                procXlength[2] += ((procsPerAxis[2]-1)==k)?(xlength+2)%procsPerAxis[2]:0;
                x1 = i*procXlength[0]; x2 = x1 + procXlength[0] + 2;
                y1 = j*procXlength[1]; y2 = y1 + procXlength[1] + 2;
                z1 = k*procXlength[2]; z2 = z1 + procXlength[2] + 2;
                snprintf(pFileName, 80, "%s.%i.%i.vts",filename,i+(j+k*procsPerAxis[1])*procsPerAxis[0],t);
                // This is stupid "../<fileName>" but what the heck
                fprintf(fp, "<Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s.%i.%i\"/>\n",x1,x2,y1,y2,z1,z2,"Debug",i+j*procsPerAxis[0]+k*procsPerAxis[0]*procsPerAxis[1],t);
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
