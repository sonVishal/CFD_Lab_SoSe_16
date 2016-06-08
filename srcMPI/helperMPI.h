
#ifndef _HELPERMPI_H_
#define _HELPERMPI_H_

#include "LBDefinitions.h"
#include <mpi/mpi.h>

void communicate(double** sendBuffer, double**readBuffer, double* collideField, const t_procData *procData);
void extract( double** sendBuffer, double* collideField, const t_procData *procData, int direction);
void swap(double** sendBuffer, double**readBuffer, const t_procData *procData, int direction);
void inject(double** readBuffer, double* collideField, const t_procData *procData, int direction);

void p_assignIndices(int *direction, int *index);

#endif
