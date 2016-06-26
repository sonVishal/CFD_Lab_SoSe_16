#ifndef _UNITTEST_H_
#define _UNITTEST_H_
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <mpi/mpi.h>

#define _TOL_ 1e-15

void storeMassVector(const t_component * const c, double ** massVector,
    const int * const xlength);

void checkMassVector(double *massVectorBefore[], double *massVectorAfter[],
    const int * const xlength, const int rank);

void computeCellMomentum(const double * const currentCell, double *momentum);

void computeGlobalMomentum(const t_component * const c,
    const int * const xlength, double * compMomentum);

void checkMomentum(const double * const momentumBefore, const double * const momentumAfter);
#endif /* end of include guard: _UNITTEST_H_ */
