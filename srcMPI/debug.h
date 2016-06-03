#include <stdio.h>

void writevtkHeaderDebug(FILE *fp, int xlength);
void writevtkPointCoordinatesDebug(FILE *fp, int xlength);

void writeVtkDebug(const double * const collideField,
    const int * const flagField, const char * filename, int xlength);

void writeCollideFieldDebug(double* collideField, char *filename, int size);
