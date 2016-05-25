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
    unsigned int t, int *xlength)
{
    // Files related variables
    char pFileName[80];
    char pTempFile[80];
    FILE *tmp = NULL;
    FILE *fp = NULL;

    // temporary variable for concatenation
    char ch;

    // Create the file with time information in the name
    sprintf(pFileName, "%s.%i.vtk", filename, t);
    sprintf(pTempFile, "temp.vtk");
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
    writevtkHeader(fp, xlength);

    // Write the point data for the domain
    writevtkPointCoordinates(fp, xlength);

    int x, y, z;            // iteration variables
    int idx;                // cell index
    double cellDensity;     // cell density

    // cell average velocity
    double cellVelocity[3];

    // Open two files and concatenate them at the end

    // Write cell velocity to the vtk file
    fprintf(fp,"\nCELL_DATA %d \n", xlength[0]*xlength[1]*xlength[2]);
    fprintf(fp, "\nVECTORS velocity float\n");

    // Write cell average density to a temporary vtk file
    fprintf(tmp, "SCALARS density float 1 \n");
    fprintf(tmp, "LOOKUP_TABLE default \n");
    for(z = 1; z <= xlength[2]; z++) {
        for(y = 1; y <= xlength[1]; y++) {
            for(x = 1; x <= xlength[0]; x++) {
                // Compute the base index for collideField
                int xyzoffset = z*(xlength[0]+2)*(xlength[1]+2) + y*(xlength[0]+2) + x;

                if (flagField[xyzoffset] == FLUID) {
                    // If it is a fluid cell then write the actual
                    // velocity and density
                    idx = Q*xyzoffset;
                    computeDensity(&collideField[idx], &cellDensity);
                    computeVelocity(&collideField[idx], &cellDensity, &cellVelocity[0]);
                } else {
                    // Else write values to signify that the cell is
                    // not a fluid cell.
                    // TODO: (VS) Decide on the values
                    cellVelocity[0] = 0.0;
                    cellVelocity[1] = 0.0;
                    cellVelocity[2] = 0.0;
                    cellDensity     = -1.0;
                }
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

void writeVtkDebug(const double * const collideField,
    const int * const flagField, const char * filename, int *xlength)
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

    int x, y, z;            // iteration variables

    fprintf(fp,"\nCELL_DATA %d \n", (xlength[0]+2)*(xlength[1]+2)*(xlength[2]+2));

    // Write cell average density to a temporary vtk file
    fprintf(fp, "SCALARS density integer 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for(z = 0; z <= xlength[2]+1; z++) {
        for(y = 0; y <= xlength[1]+1; y++) {
            for(x = 0; x <= xlength[0]+1; x++) {
                // Compute the base index for collideField
                int xyzoffset = z*(xlength[0]+2)*(xlength[1]+2) + y*(xlength[0]+2) + x;
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

// Header for the VTK files
void writevtkHeader(FILE *fp, int *xlength)
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
    fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength[0]+1, xlength[1]+1, xlength[2]+1);
    fprintf(fp,"POINTS %i integer\n", (xlength[0]+1)*(xlength[1]+1)*(xlength[2]+1));
    fprintf(fp,"\n");
}

void writevtkHeaderDebug(FILE *fp, int *xlength)
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
    fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength[0]+3, xlength[1]+3, xlength[2]+3);
    fprintf(fp,"POINTS %i integer\n", (xlength[0]+3)*(xlength[1]+3)*(xlength[2]+3));
    fprintf(fp,"\n");
}

void writevtkPointCoordinates(FILE *fp, int *xlength) {
    int x, y, z;

    // We have xlength + 1 points for xlength cells in each direction
    for(z = 0; z <= xlength[2]; z++) {
        for(y = 0; y <= xlength[1]; y++) {
            for(x = 0; x <= xlength[0]; x++) {
                fprintf(fp, "%d %d %d\n", x, y, z);
            }
        }
    }
}

void writevtkPointCoordinatesDebug(FILE *fp, int *xlength) {
    int x, y, z;

    // We have xlength + 3 points for xlength+2 cells in each direction
    for(z = 0; z <= xlength[2]+2; z++) {
        for(y = 0; y <= xlength[1]+2; y++) {
            for(x = 0; x <= xlength[0]+2; x++) {
                fprintf(fp, "%d %d %d\n", x, y, z);
            }
        }
    }
}
