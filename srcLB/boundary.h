#include "LBDefinitions.h"
#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

/* handles the boundaries in our simulation setup */
void treatBoundary(double *collideField, const t_flagField * const flagField,
	const t_boundPara * const boundPara, const int * const xlength);

void p_noSlip(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	int const * const totalSize);

void p_movingWall(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize);

void p_freeSlip(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	int const * const totalSize);

void p_outflow(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize);

void p_inflow(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara);

void p_pressureIn(double* collideField, t_flagField const * const flagField,
	int const * const point, int const * const xlength,
	const t_boundPara * const boundPara, int const * const totalSize);

void p_assignIndices(const short int * const flag, int * index, int * mirrorIndex);
#endif
