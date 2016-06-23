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
void writeVtsOutput(const t_component * const c, const int numComp,
    const int * const flagField, const char * filename, const unsigned int t,
    const int xlen, const t_procData * const procData, const int * const procsPerAxis);

void writevtsPointCoordinates(FILE * fp, const int xlen,
    const int * const xlength, const int * const myPos, const int * const procsPerAxis);

void p_writeCombinedPVTSFile(const int numComp, const char * filename,
    const unsigned int t, const int xlength, const int * const procsPerAxis);
#endif
