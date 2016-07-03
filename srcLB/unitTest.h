#ifndef _UNITTEST_H_
#define _UNITTEST_H_
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <mpi/mpi.h>
#include "stdio.h"
#include "computeCellValues.h"

#define _TOL_ 1e-15

void storeMassVector(const t_component * const c, const int xlength, double **massVector);

void checkMassVector(const int xlength);

void computeCellMomentum(const double * const currentCell, double *momentum);

void initializeUnitTest(const int ts);

void computeGlobalMomentum(const t_component * const c, const int xlength, double *momentum);

void beforeCollision(const t_component * const c, int xLength);

void afterCollision(const t_component * const c, int xLength);

void checkMomentum();

void freeUnitTest();

#endif /* end of include guard: _UNITTEST_H_ */
