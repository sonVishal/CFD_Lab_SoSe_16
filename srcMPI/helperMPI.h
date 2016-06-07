
#include "LBDefinitions.h"
#include <mpi/mpi.h>

void communicate(double** sendBuffer, double**readBuffer, double* collideField, int* xlength, int* bufferSize);
void extract( double** sendBuffer, double* collideField, int* xlength, int* bufferSize, int direction);
void swap(double** sendBuffer, double**readBuffer, int* bufferSize, int direction);
void inject(double** readBuffer, double* collideField, int* xlength, int* bufferSize, int direction);
