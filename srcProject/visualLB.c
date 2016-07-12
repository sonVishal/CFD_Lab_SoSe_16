#include "visualLB.h"

// Function to write VTS output files for visualization
void writeVtsOutput(const t_component * const c, const char * filename,
    unsigned int t, int xlen, const t_procData * const procData,
    const int * const procsPerAxis) {
    // Files related variables
    char pFileName[80];
    FILE *fp = NULL;

    // Create the file with time information in the name
    sprintf(pFileName, "%s.%i.%i.vts", filename,procData->rank,t);
    fp  = fopen(pFileName, "w");

    int myPos[3] = {0,0,0};
    // printf("Rank %d\n",procData->rank);
    p_rankToPos(procsPerAxis, procData->rank, myPos);

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
    writevtsPointCoordinates(fp,xlen,procData->xLength,myPos,procsPerAxis);

    fprintf(fp,"<CellData>\n");

    for (int i = 0; i < numComp; i++) {
        fprintf(fp,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity_C%d\">\n",i);

        int x, y, z;            // iteration variables
        int idx;                // cell index
        double cellDensity;     // cell density

        // cell average velocity
        double cellVelocity[3] = {0.0,0.0,0.0};

        // Temporary variables for (xlength+2)^2
        int const xylen = (procData->xLength[0]+2)*(procData->xLength[1]+2);

        // Temporary variables for z and y offsets
        int zOffset, yzOffset;

        // Write cell average velcity the file
        for(z = 1; z <= procData->xLength[2]; z++) {
            zOffset = z*xylen;
            for(y = 1; y <= procData->xLength[1]; y++) {
                yzOffset = zOffset + y*(procData->xLength[0]+2);
                for(x = 1; x <= procData->xLength[0]; x++) {
                    // Compute the base index for field_new
                    idx = Q*(yzOffset + x);

                    computeNumDensity(&c[i].field_new[idx], &cellDensity);
                    computeVelocity(&c[i].field_new[idx], &cellDensity, cellVelocity);

                    // Write cell average velocities
                    fprintf(fp, "%f %f %f\n", cellVelocity[0],
                    cellVelocity[1], cellVelocity[2]);
                }
            }
        }
        fprintf(fp,"</DataArray>\n");
        fprintf(fp,"<DataArray type=\"Float32\" NumberOfComponents=\"1\" Name=\"Density_C%d\">\n",i);
        // Write cell average denisties the file
        for(z = 1; z <= procData->xLength[2]; z++) {
            zOffset = z*xylen;
            for(y = 1; y <= procData->xLength[1]; y++) {
                yzOffset = zOffset + y*(procData->xLength[0]+2);
                for(x = 1; x <= procData->xLength[0]; x++) {
                    // Compute the base index for field_new
                    idx = Q*(yzOffset + x);

                    computeNumDensity(&c[i].field_new[idx], &cellDensity);
                    // Write cell average velocities
                    fprintf(fp, "%f\n", cellDensity);
                }
            }
        }
        fprintf(fp,"</DataArray>\n");
    }

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

// Function to write coordinates of the points
void writevtsPointCoordinates(FILE * fp, const int xlen, const int * const xlength,
    const int * const myPos, const int * const procsPerAxis) {
    int x, y, z;
    int baseLength[3] = {xlen/procsPerAxis[0], xlen/procsPerAxis[1], xlen/procsPerAxis[2]};
    // If the proc is at the end of some axis then add the remaining length
    int myLen[3] = {0,0,0};
    myLen[0] = (myPos[0] == procsPerAxis[0]-1)?xlen-(procsPerAxis[0]-1)*baseLength[0]:baseLength[0];
    myLen[1] = (myPos[1] == procsPerAxis[1]-1)?xlen-(procsPerAxis[1]-1)*baseLength[1]:baseLength[1];
    myLen[2] = (myPos[2] == procsPerAxis[2]-1)?xlen-(procsPerAxis[2]-1)*baseLength[2]:baseLength[2];
    unsigned int x1 = myPos[0]*baseLength[0];
    unsigned int x2 = x1 + myLen[0];
    unsigned int y1 = myPos[1]*baseLength[1];
    unsigned int y2 = y1 + myLen[1];
    unsigned int z1 = myPos[2]*baseLength[2];
    unsigned int z2 = z1 + myLen[2];
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

// Function to combine the VTS files into one
void p_writeCombinedPVTSFile(const char * const filename, const unsigned int t,
    const int xlen, const int * const procsPerAxis) {
    // Files related variables
    char pFileName[80];
    FILE *fp = NULL;

    // Create the file with time information in the name
    sprintf(pFileName, "%s.%i.pvts", "pv_files/project_combined",t);
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
    fprintf(fp, "<PStructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\" GhostLevel=\"1\">\n",xlen,xlen,xlen);
    fprintf(fp, "<PPoints>\n");
    fprintf(fp, "%s\n","<PDataArray NumberOfComponents=\"3\" type=\"UInt32\" />");
    fprintf(fp, "</PPoints>\n");
    fprintf(fp, "<PCellData>\n");
    for (int i = 0; i < numComp; i++) {
        fprintf(fp, "<PDataArray type=\"Float32\" NumberOfComponents=\"3\" Name=\"Velocity_C%d\"/>\n",i);
        fprintf(fp, "<PDataArray type=\"Float32\" NumberOfComponents=\"1\" Name=\"Density_C%d\"/>\n",i);
    }
    fprintf(fp, "</PCellData>\n");

    int baseLength[3] = {xlen/procsPerAxis[0], xlen/procsPerAxis[1], xlen/procsPerAxis[2]};

    int myLen[3] = {0,0,0};

    for (int k = 0; k < procsPerAxis[2]; k++) {
        for (int j = 0; j < procsPerAxis[1]; j++) {
            for (int i = 0; i < procsPerAxis[0]; i++) {
                myLen[0] = (i == procsPerAxis[0]-1)?xlen-(procsPerAxis[0]-1)*baseLength[0]:baseLength[0];
                myLen[1] = (j == procsPerAxis[1]-1)?xlen-(procsPerAxis[1]-1)*baseLength[1]:baseLength[1];
                myLen[2] = (k == procsPerAxis[2]-1)?xlen-(procsPerAxis[2]-1)*baseLength[2]:baseLength[2];
                unsigned int x1 = i*baseLength[0];
                unsigned int x2 = x1 + myLen[0];
                unsigned int y1 = j*baseLength[1];
                unsigned int y2 = y1 + myLen[1];
                unsigned int z1 = k*baseLength[2];
                unsigned int z2 = z1 + myLen[2];
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
