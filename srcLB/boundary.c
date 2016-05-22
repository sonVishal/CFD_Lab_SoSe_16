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

void p_handleInflow(){
	/* TODO: See page 12 for description and equation (2.2) */
	ERROR("TODO");
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

	for(int z = 0; z <= xlength[0]+1; ++z){
		for(int y = 0; y <= xlength[1]+1; ++y){
			for(int x = 0; x <= xlength[2]+1; ++x){
				//TODO: (DL) make function for index (therefore not optimized)
				xyz_offset = (z*xlength[2]*xlength[2] + y*xlength[1] + x);
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
					case INFLOW:
						p_handleInflow();
						break;
					case OUTFLOW:
						p_handleOutflow();
						break;
					case PRESSURE_IN:
						p_handleInPressure();
						break;
					case OBSTACLE:
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
