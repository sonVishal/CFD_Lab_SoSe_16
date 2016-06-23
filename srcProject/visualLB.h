#ifndef _VISUALLB_H_
#define _VISUALLB_H_
#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "helper.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by 'filename' and timestep 't'. You can re-use parts of the code
 *  from visual.c (VTK output for Navier-Stokes solver) and modify it for 3D datasets.
 */
void writeVtsOutput(const t_component * const c, int numComp, const int * const flagField, const char * filename, unsigned int t, int xlen, t_procData procData, int *procsPerAxis);

void writevtsPointCoordinates(FILE *fp, int xlen, int *xlength, int *myPos, int *procsPerAxis);

void p_writeCombinedPVTSFile(const int numComp, const char * filename, unsigned int t, int xlength, int *procsPerAxis);
#endif
