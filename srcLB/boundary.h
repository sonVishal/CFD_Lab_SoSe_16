#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "LBDefinitions.h"

/* handles the boundaries in our simulation setup */
void treatBoundary(double *collideField, const int * const flagField,
	const t_boundPara * const boundPara, const int * const xlength);

void p_noSlip(double* collideField, int const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const normal,
	int const * const totalSize);

void p_movingWall(double* collideField, int const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const normal,
	int const * const totalSize);

void p_freeSlip(double* collideField, int const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const normal,
	int const * const totalSize);

void p_outflow(double* collideField, int const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const normal,
	int const * const totalSize);

void p_pressureIn(double* collideField, int const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const normal,
	int const * const totalSize);

void p_assignIndices(const int * const normal, int * index, int * mirrorIndex);
#endif
