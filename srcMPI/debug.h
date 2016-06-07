#include <stdio.h>
#include "LBDefinitions.h"

void writevtkHeaderDebug(FILE *fp, int xlength);
void writevtkPointCoordinatesDebug(FILE *fp, int xlength);

void writeVtkDebug(const double * const collideField,
    const int * const flagField, const char * filename, int xlength);

void writeCollideFieldDebug(char *filename, double* collideField, int size);

void checkCollideFieldDebug(char *referenceFile, double *currentCollideField, int size);

void printProcData(t_procData procData);

void printProcDataPos(t_procData procData, int *pos);
