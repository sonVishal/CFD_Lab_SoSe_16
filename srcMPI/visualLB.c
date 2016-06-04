#include "visualLB.h"
#include "computeCellValues.h"
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
    unsigned int t, int xlength, int rank, int number_of_ranks)
{
    // Files related variables
    char pFileName[80];
    char pTempFile[80];
    FILE *tmp = NULL;
    FILE *fp = NULL;

    // temporary variable for concatenation
    char ch;

    // Create the file with time information in the name
    sprintf(pFileName, "%s.%i.%i.vtk", filename, t, rank);
    sprintf(pTempFile, "temp.%i.vtk", rank);
    fp  = fopen(pFileName, "w");
    tmp = fopen(pTempFile, "w");

    // Check if files were opened or not
    if(fp == NULL)
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", pFileName);
        ERROR(szBuff);
        return;
    }
    if(tmp == NULL)
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", pTempFile);
        ERROR(szBuff);
        return;
    }

    // Write header for the VTK file
    writevtkHeader(fp,xlength);

    // Write the point data for the domain
    writevtkPointCoordinates(fp,xlength);

    int x, y, z;            // iteration variables
    int idx;                // cell index
    double cellDensity;     // cell density

    // cell average velocity
    double cellVelocity[3] = {0.0,0.0,0.0};

    // Temporary variables for (xlength+2)^2
    int const xlen2 = (xlength+2)*(xlength+2);

    // Temporary variables for z and y offsets
    int zOffset, yzOffset;

    // Open two files and concatenate them at the end

    // Write cell velocity to the vtk file
    fprintf(fp,"\nCELL_DATA %d \n", xlength*xlength*xlength);
    fprintf(fp, "\nVECTORS velocity float\n");

    // Write cell average density to a temporary vtk file
    fprintf(tmp, "SCALARS density float 1 \n");
    fprintf(tmp, "LOOKUP_TABLE default \n");
    for(z = 1; z <= xlength; z++) {
        zOffset = z*xlen2;
        for(y = 1; y <= xlength; y++) {
            yzOffset = zOffset + y*(xlength+2);
            for(x = 1; x <= xlength; x++) {
                // Compute the base index for collideField
                idx = Q*(yzOffset + x);

                computeDensity(&collideField[idx], &cellDensity);
                computeVelocity(&collideField[idx], &cellDensity, &cellVelocity[0]);

                // Write cell average velocities
                fprintf(fp, "%f %f %f\n", cellVelocity[0],
                cellVelocity[1], cellVelocity[2]);

                // Write the cell density
                fprintf(tmp, "%f\n", cellDensity);
            }
        }
    }

    // Close file
    if(fclose(tmp))
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", pTempFile);
        ERROR(szBuff);
    }

    // Open for reading
    tmp = fopen(pTempFile, "r");
    if(tmp == NULL)
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", pTempFile);
        ERROR(szBuff);
        return;
    }

    // Concatenate the two files
    while((ch = fgetc(tmp)) != EOF)
        fputc(ch,fp);

    // Close files
    if(fclose(fp))
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", pFileName);
        ERROR(szBuff);
    }
    if(fclose(tmp))
    {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", pTempFile);
        ERROR(szBuff);
    }

    // Delete the temporary file
    if (remove(pTempFile)) {
        char szBuff[80];
        sprintf(szBuff, "Failed to delete %s", pTempFile);
        ERROR(szBuff);
    }
}

// Header for the VTK files
void writevtkHeader(FILE *fp, int xlength)
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
    fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength+1, xlength+1, xlength+1);
    fprintf(fp,"POINTS %i integer\n", (xlength+1)*(xlength+1)*(xlength+1));
    fprintf(fp,"\n");
}

void writevtkPointCoordinates(FILE *fp, int xlength) {
    int x, y, z;

    // We have xlength + 1 points for xlength cells in each direction
    for(z = 0; z <= xlength; z++) {
        for(y = 0; y <= xlength; y++) {
            for(x = 0; x <= xlength; x++) {
                fprintf(fp, "%d %d %d\n", x, y, z);
            }
        }
    }
}
