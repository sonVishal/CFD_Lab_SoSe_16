#ifndef _DEBUG_H_
#define _DEBUG_H_

#if _POSIX_C_SOURCE < 2
    #define _POSIX_C_SOURCE 3
#endif

#include "helper.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include "LBDefinitions.h"

void writeVtkOutputDebug(const double * const collideField, const int * const flagField, const char * filename, unsigned int t, t_procData procData, int *procsPerAxis);

void writevtkPointCoordinatesDebug(FILE *fp, int *xlength, int *myPos);

void p_writeCombinedPVTSFileDebug(const char * filename, unsigned int t, int xlength, int *procsPerAxis);

void writeCollideFieldDebug(char *filename, double* collideField, int size);

void checkCollideFieldDebug(char *referenceFile, double *currentCollideField, int size);

void printProcData(t_procData procData);

void printProcDataPos(t_procData procData, int *pos);

void convertEnumWallToString(const int wall, char *wallName);

void printWallEnum(const int wall);

void convertEnumCellToString(const int wall, char *wallName);

void printCellEnum(const int wall);

void debug_setBufferValues(double **sendBuffer, double **readBuffer, t_procData procData);

int parse_output(void);

#endif
