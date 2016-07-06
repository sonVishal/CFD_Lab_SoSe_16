#ifndef _UNITTEST_H_
#define _UNITTEST_H_
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include "parallel.h"
#include <mpi/mpi.h>

#define _TOL_ 1e-15

void storeMassVector(const t_component * const c, const int * const xlength, double **massVector);

void checkMassVector(const int * const xlength, const int rank);

void computePostCollisionDistributions(const t_component * const c, t_procData const * const procData);

void doStreaming(const t_component * const c, t_procData const * const procData);

void steamCollideUnitTest(const t_component * const c, t_procData const * const procData);

void treatPostCollisionBoundary(const t_component * const c, t_procData const * const procData, int direction);

void computeCellMomentum(const double * const currentCell, double *momentum);

void initializeUnitTest(const int ts);

void computeGlobalMomentum(const t_component * const c, const int * const xlength, double *momentum);

void beforeCollision(const t_component * const c, t_procData const * const procData);

void afterCollision(const t_component * const c, t_procData const * const procData);

void checkMomentum();

void freeUnitTest();

#endif /* end of include guard: _UNITTEST_H_ */
