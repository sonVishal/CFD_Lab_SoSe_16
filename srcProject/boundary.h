//#ifndef _BOUNDARY_H_
//#define _BOUNDARY_H_
//#include "parallel.h"
//#include "helper.h"
//#include "LBDefinitions.h"
//#include "computeCellValues.h"
//#include <stdio.h>
//#include <mpi/mpi.h>


//typedef struct{
    //int x, y, z;
//} t_iterParaEdge;

//[>Wrapper around tratBoundary to handle every component<]
//void treatComponentBoundary(t_component *c, int const * const flagField, const t_procData * const procData, double **sendBuffer, double **readBuffer, int densityFlag);

//[> handles the boundaries in our simulation setup <]
//void treatBoundary(int const * const flagField, double *collideField, const t_procData * const procData, double **sendBuffer, double **readBuffer, const int densityFlag);

////TODO (DL) make comments for function
//void treatPeriodicEdgeNoComm(double *collideField, const t_procData * const procData, const int edge1, const int edge2, const int densityFlag);
//void treatPeriodicWallNoComm(double *collideField, const int wall1, const int wall2, t_procData const*const procData, const int densityFlag);
//void treatPeriodicEdge(double *collideField, double *const sendBuffer, double *const readBuffer, const t_procData * const procData, const int procEdge, const int opponentEdge, const int densityFlag);
//void treatPeriodicWall(int const * const flagField,
    //double *const collideField, double *const sendBuffer, double *const readBuffer,
    //const t_procData * const procData, const int procWall, const int opponentWall, const int densityFlag);
//int p_assignSharedEdgeIndex(const int edge);
//void p_setEdgeIterParameters(t_iterParaEdge *const iterPara, t_procData const*const procData, const int edge, const int injectFlag);

////TODO: (DL) I get compiler errors for this because of t_iterPara, even though parallel.h is included it does not know it...
//// void p_setBoundaryIterParameters(t_iterPara *const iterPara, t_procData const*const procData, const int direction);


//int extractInjectEdge(double buffer[], double * const collideField, t_iterParaEdge const * const iterPara, t_procData const * const procData, const int index, const int injectFlag, const int densityFlag);

//#endif


#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "LBDefinitions.h"
#include "helper.h"
#include "computeCellValues.h"
#include <stdio.h>

/** handles the boundaries in our simulation setup */
void treatBoundary(t_component *collideField, int xlength);
int computeCellOffset(const int outer, const int inner, const int fixed, const int dir, const int xlength);
void treatWallPeriodic(t_component * collideField, int direction, int xlength);
#endif
