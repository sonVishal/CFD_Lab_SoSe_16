#include "debug.h"
#include "helper.h"
#include <stdio.h>
#include "LBDefinitions.h"

void writeVtkDebug(const double * const collideField,
    const int * const flagField, const char * filename, int xlength)
{
    // Files related variables
    char pFileName[80];
    FILE *fp = NULL;

    // Create the file with time information in the name
    sprintf(pFileName, "%s.vtk", filename);
    fp  = fopen(pFileName, "w");

    // Check if files were opened or not
    if(fp == NULL)
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", pFileName);
        ERROR(szBuff);
        return;
    }

    // Write header for the VTK file
    writevtkHeaderDebug(fp, xlength);

    // Write the point data for the domain
    writevtkPointCoordinatesDebug(fp, xlength);

    // iteration variables
    int x, y, z;

    fprintf(fp,"\nCELL_DATA %d \n", (xlength+2)*(xlength+2)*(xlength+2));

    // Write cell average density to a temporary vtk file
    fprintf(fp, "SCALARS boundaryType integer 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for(z = 0; z <= xlength+1; z++) {
        for(y = 0; y <= xlength+1; y++) {
            for(x = 0; x <= xlength+1; x++) {
                // Compute the base index for collideField
                int xyzoffset = z*(xlength+2)*(xlength+2) + y*(xlength+2) + x;
                fprintf(fp, "%d\n", flagField[xyzoffset]);
            }
        }
    }
    // Close files
    if(fclose(fp))
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", pFileName);
        ERROR(szBuff);
    }
}

void writevtkHeaderDebug(FILE *fp, int xlength)
{
    if(fp == NULL)
    {
        char szBuff[80];
        sprintf(szBuff, "Null pointer in write_vtkHeader");
        ERROR(szBuff);
        return;
    }

    fprintf(fp,"# vtk DataFile Version 2.0\n");
    fprintf(fp,"generated by CFD-lab course output (written by Vishal Sontakke) \n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"\n");
    fprintf(fp,"DATASET STRUCTURED_GRID\n");
    fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength+3, xlength+3, xlength+3);
    fprintf(fp,"POINTS %i integer\n", (xlength+3)*(xlength+3)*(xlength+3));
    fprintf(fp,"\n");
}

void writevtkPointCoordinatesDebug(FILE *fp, int xlength) {
    int x, y, z;

    // We have xlength + 3 points for xlength+2 cells in each direction
    for(z = 0; z <= xlength+2; z++) {
        for(y = 0; y <= xlength+2; y++) {
            for(x = 0; x <= xlength+2; x++) {
                fprintf(fp, "%d %d %d\n", x, y, z);
            }
        }
    }
}


void writeCollideFieldDebug(char *filename, double* collideField, int size){

    FILE *fp = NULL;

    char filenameExt[100];
    snprintf(filenameExt, 100, "%s.array", filename);

    fp  = fopen(filenameExt, "w");

    if(fp == NULL)
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", filenameExt);
        ERROR(szBuff);
        return;
    }

    int linebreak = Q;
    int counter = 0;
    int keepWriting = 1;

    while(keepWriting){
    	for(int l = 0; l < linebreak && keepWriting; ++l){
    		fprintf(fp,"%.16lf ", collideField[counter]);
    		counter++;

    		if(counter >= size-1)
    			keepWriting = 0;
    	}
    	fprintf(fp, "\n");
    }

    if(fclose(fp))
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", filenameExt);
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