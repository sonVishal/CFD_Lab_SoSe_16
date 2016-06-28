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
#include "computeCellValues.h"

void writeVtsOutputDebug(const t_component * const c, const int * const flagField,
    const char * filename, unsigned int t, int xlen,
    const t_procData * const procData, const int * const procsPerAxis);

void writevtsPointCoordinatesDebug(FILE * fp, const int xlen, const int * const xlength,
    const int * const myPos, const int * const procsPerAxis);

void p_writeCombinedPVTSFileDebug(const char * const filename, const unsigned int t,
    const int xlen, const int * const procsPerAxis);

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
