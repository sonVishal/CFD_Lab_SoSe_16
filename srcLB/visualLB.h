#include <stdio.h>
#include "LBDefinitions.h"
#include <unistd.h>

#ifndef _VISUALLB_H_
#define _VISUALLB_H_

/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by 'filename' and timestep 't'. You can re-use parts of the code
 *  from visual.c (VTK output for Navier-Stokes solver) and modify it for 3D datasets.
 */
void writeVtkOutput(const double * const collideField, const int * const flagField, const char * filename, unsigned int t, int *xlength);
void writeVtkDebug(const double * const collideField,
    const int * const flagField, const char * filename, int *xlength);

void writevtkPointCoordinates(FILE *fp, int *xlength);
void writevtkPointCoordinatesDebug(FILE *fp, int *xlength);

void writevtkHeader(FILE *fp, int *xlength);

#endif
