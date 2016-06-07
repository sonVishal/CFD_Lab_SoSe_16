#include "visualLB.h"
#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "helper.h"
#include <stdlib.h>

// Function to write VTK output files for visualization
// Inputs
// collideField - Probability distribution function
//                  Length = Q*(xlength+2)^D
// flagField    - Geometry info about whether a cell is
//                  FLUID         = 0
//                  NO_SLIP       = 1
//                  MOVING_WALL   = 2
//                  Length        = (xlength+2)^D
// filename     - Name of the file to which output is to be written
// t            - Time at which output is to be stored
// xlength      - Number of cells in one direction

void writeVtkOutput(const double * const collideField,
    const int * const flagField, const char * filename,
    unsigned int t, t_procData procData, int *procsPerAxis)
{
    // Files related variables
    char pFileName[80];
    FILE *fp = NULL;

    // Create the file with time information in the name
    sprintf(pFileName, "%s.%i.%i.vts", filename,procData.rank,t);
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
    writevtkPointCoordinates(fp,procData.xLength,myPos);

    fprintf(fp,"<CellData>\n");
    fprintf(fp,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\">\n");

    int x, y, z;            // iteration variables
    int idx;                // cell index
    double cellDensity;     // cell density

    // cell average velocity
    double cellVelocity[3] = {0.0,0.0,0.0};

    // Temporary variables for (xlength+2)^2
    int const xylen = (procData.xLength[0]+2)*(procData.xLength[1]+2);

    // Temporary variables for z and y offsets
    int zOffset, yzOffset;

    // Write cell average velcity the file
    for(z = 1; z <= procData.xLength[2]; z++) {
        zOffset = z*xylen;
        for(y = 1; y <= procData.xLength[1]; y++) {
            yzOffset = zOffset + y*(procData.xLength[0]+2);
            for(x = 1; x <= procData.xLength[0]; x++) {
                // Compute the base index for collideField
                idx = Q*(yzOffset + x);

                computeDensity(&collideField[idx], &cellDensity);
                computeVelocity(&collideField[idx], &cellDensity, &cellVelocity[0]);

                // Write cell average velocities
                fprintf(fp, "%f %f %f\n", cellVelocity[0],
                cellVelocity[1], cellVelocity[2]);
            }
        }
    }
    fprintf(fp,"</DataArray>\n");
    fprintf(fp,"<DataArray type=\"Float32\" NumberOfComponents=\"1\" Name=\"Density\">\n");
    // Write cell average denisties the file
    for(z = 1; z <= procData.xLength[2]; z++) {
        zOffset = z*xylen;
        for(y = 1; y <= procData.xLength[1]; y++) {
            yzOffset = zOffset + y*(procData.xLength[0]+2);
            for(x = 1; x <= procData.xLength[0]; x++) {
                // Compute the base index for collideField
                idx = Q*(yzOffset + x);

                computeDensity(&collideField[idx], &cellDensity);
                // Write cell average velocities
                fprintf(fp, "%f\n", cellDensity);
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

void writevtkPointCoordinates(FILE *fp, int *xlength, int *myPos) {
    int x, y, z;
    // printf("Position = (%d,%d,%d)\n",myPos[0],myPos[1],myPos[2]);
    unsigned int x1 = myPos[0]*(xlength[0]);
    unsigned int x2 = x1 + xlength[0];
    unsigned int y1 = myPos[1]*(xlength[1]);
    unsigned int y2 = y1 + xlength[1];
    unsigned int z1 = myPos[2]*(xlength[2]);
    unsigned int z2 = z1 + xlength[2];
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

void p_writeCombinedPVTSFile(const char * filename, unsigned int t, int xlength, int *procsPerAxis) {
    // Files related variables
    char pFileName[80];
    FILE *fp = NULL;

    // Create the file with time information in the name
    sprintf(pFileName, "%s.%i.pvts", "pv_files/worksheet4",t);
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
    fprintf(fp, "<PStructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\" GhostLevel=\"1\">\n",xlength,xlength,xlength);
    fprintf(fp, "<PPoints>\n");
    fprintf(fp, "%s\n","<PDataArray NumberOfComponents=\"3\" type=\"UInt32\" />");
    fprintf(fp, "</PPoints>\n");
    fprintf(fp, "<PCellData>\n");
    fprintf(fp, "<PDataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity\"/>\n");
    fprintf(fp, "<PDataArray type=\"Float32\" NumberOfComponents=\"1\" Name=\"Density\"/>\n");
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
                procXlength[0] += ((procsPerAxis[0]-1)==i)?xlength%procsPerAxis[0]:0;
                procXlength[1] += ((procsPerAxis[1]-1)==j)?xlength%procsPerAxis[1]:0;
                procXlength[2] += ((procsPerAxis[2]-1)==k)?xlength%procsPerAxis[2]:0;
                x1 = i*procXlength[0]; x2 = x1 + procXlength[0];
                y1 = j*procXlength[1]; y2 = y1 + procXlength[1];
                z1 = k*procXlength[2]; z2 = z1 + procXlength[2];
                snprintf(pFileName, 80, "%s.%i.%i.vts",filename,i+(j+k*procsPerAxis[1])*procsPerAxis[0],t);
                // This is stupid "../<fileName>" but what the heck
                fprintf(fp, "<Piece Extent=\"%d %d %d %d %d %d\" Source=\"../%s\"/>\n",x1,x2,y1,y2,z1,z2,pFileName);
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
