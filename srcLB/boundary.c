#include <stdio.h>
#include "boundary.h"
#include "helper.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void p_handleNoSlip(double *collideField, int const * const flagField, const int current_cell_index,
		int const * const point, int const * const xlength){
	int n_cell_index;

	for(int i = 0; i < Q; i++){
		int n_point[3] = {point[0] + LATTICEVELOCITIES[i][0],
				point[1] + LATTICEVELOCITIES[i][1],
				point[2] + LATTICEVELOCITIES[i][2]};
		n_cell_index = Q*(n_point[2]*xlength[2]*xlength[2] + n_point[1]*xlength[1] + n_point[0]);

		if(flagField[n_cell_index] == FLUID){ /* TODO: (DL) check also for valid index */
			collideField[current_cell_index + i] = collideField[n_cell_index + (Q-i-1)];
		}
	}
}

void p_handleMovingWall(double *collideField, const int currentCellIndex,
		int const * const point, int const * const xlength, double const * const wallVelocity){

	int n_cell_index;
	double density;

	for(int i = 0; i < Q; i++){
		int c[3] = {LATTICEVELOCITIES[i][0], LATTICEVELOCITIES[i][1], LATTICEVELOCITIES[i][2]};
		int n_point[3] = {point[0] + c[0], point[1] + c[1], point[2] + c[2]};
		n_cell_index = Q*(n_point[2]*xlength[2]*xlength[2] + n_point[1]*xlength[1] + n_point[0]);

		if(n_cell_index == FLUID){ /* TODO: (DL) check also for valid index */
			double dot_uwall_c = wallVelocity[0]*c[0]+wallVelocity[1]*c[1]+wallVelocity[2]*c[2];
			double weight = LATTICEWEIGHTS[i];
			computeDensity(&collideField[n_cell_index], &density);

			collideField[currentCellIndex + i] = collideField[n_cell_index + (Q-i-1)] +
					2 * weight * dot_uwall_c / (C_S*C_S);
		}
	}
}


void p_handleFreeSlip(){
    //TODO: (TKS) Need to check if a boundary cell is already updated by e.g. no-slip
    //      and then add the contribution from free-slip
	ERROR("TODO");
	/*
	 * The free-slip boundaries are closely related to the no-slip boundary,
	 * as it can be seen in Fig. 2.1. Note that the implementation of these
	 * boundaries requires more effort and communication, since the particle
	 * distribution functions are not copied into the ones opposing them in
	 * the boundaries but they are mirrored into the boundary cell. By doing
	 * so the velocity tangential to the boundary is preserved. Care has to be
	 * taken at corners, since there the “mirrored” particle distributions will
	 * pass the boundary.
	 */
}

//TODO: (TKS) Does only need to b set once since no new distributions are streamed in here.
void p_handleInflow(int x, int y, int z, int *xlength, t_boundPara *boundPara,
                    double *collideField, const int currentCellIndex){
	/* TODO: See page 12 for description and equation (2.2) */
	ERROR("TODO");
	for(int i = 0; i < Q; i++){
        int flag;
        double feq[19];

        if(x == 0)
            flag = YZ_BOTTOM;
        else if(x == xlength[0]+2)
            flag = YZ_TOP;
        else if(y == 0)
            flag = XZ_BACK;
        else if(x == xlength[1]+2)
            flag = XZ_FRONT;
        else if(z == 0)
            flag = XY_LEFT;
        else if(x == xlength[2]+2)
            flag = XY_RIGHT;
        else{
            flag = -1;
            ERROR("Inflow not a boundary (Aborting)");
        }

        //TODO: (TKS) Add case for obstacle inside the domain.

        computeFeq(&boundPara[flag].rhoRef, boundPara[flag].wallVelocity , feq);
	    for(int i = 0; i < Q; i++){
            collideField[currentCellIndex + i] = feq[i];
        }
    }

}

void p_handleOutflow(){
	/* TODO: See page 11 for description and equation (2.1) */

	ERROR("TODO");
}

void p_handleInPressure(){
	/* See page 12 for description and equation (2.2) */
	ERROR("TODO");
}

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int *xlength){

	int xyz_offset, flag ,currentCellIndex;
	int const xlen2 = xlength[0]+2;
	int const ylen2 = xlength[1]+2;
	int const zlen2 = xlength[1]+2;

	for(int z = 0; z < zlen2; ++z){
		int offset1 = z*xlen2*ylen2;
		for(int y = 0; y < ylen2; ++y){
			int offset2 = offset1 + y*xlen2;
			for(int x = 0; x < xlen2; ++x){
				xyz_offset = offset2 + x;

				flag = flagField[xyz_offset];
				if (flag != FLUID) {
					currentCellIndex = Q*xyz_offset;
					int point[3] = {x,y,z};

					switch(flag){
					case NO_SLIP:
						p_handleNoSlip(collideField, flagField, currentCellIndex, point, xlength);
						break;
					case MOVING:
						p_handleMovingWall(collideField, currentCellIndex, point, xlength, wallVelocity);
						break;
					case FREE_SLIP:
						p_handleFreeSlip();
						break;
					case INFLOW: // Set in initialize fields.
						//p_handleInflow();
						break;
					case OUTFLOW:
						p_handleOutflow();
						break;
					case PRESSURE_IN:
						p_handleInPressure();
						break;
					case OBSTACLE:
						/*TODO: (DL) evaluate first if the cell is even 'connected' to a FLUID cell,
						 * if not: nothing needs to be done.
						 */
						p_handleNoSlip(collideField, flagField, currentCellIndex, point, xlength);
						break;
					default:
						ERROR("TODO");
						break;
					}
				}
			}
		}
	}
}
