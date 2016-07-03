
#ifndef _STREAMING_H_
#define _STREAMING_H_
#include "LBDefinitions.h"

void streamCollide(t_component* c, int xlength, double* feq, int* flagField);
void updateFeq(const int *xlength, const double*rho, double *velocity[3], double*feq);

#endif
