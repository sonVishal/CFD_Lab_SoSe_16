#include "uvp.h"
#include <math.h>
#include <float.h>


void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
){
	/* TODO: implement */
}

void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
){
	/* TODO: implement */
}



double p_find_max(double **matrix, int imax, int jmax){
	int i,j;

	/* equals: minimum (negative) double value; note the "-" (!!)
	 * https://stackoverflow.com/questions/1153548/minimum-double-value-in-c-c*/
	double found_max = -DBL_MAX;


	/*TODO: (DL) this is most likely wrong, check out again how the matrix is structured */
	for(i = 0; i<imax; ++i){
		for(j = 0; j < jmax; ++j){
			if(found_max < matrix[i][j]){
				found_max = matrix[i][j];
			}
		}
	}
	return found_max;
}

void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
){
	double restriction; /* holds the currently tested restriction */
	double min;
	restriction = Re/(2 * (1/(dx*dx) + 1/(dy*dy)));
	min = restriction;

	restriction = dx / p_find_max(U, imax, jmax);
	if(min > restriction){
		min = restriction;
	}

	restriction = dy / p_find_max(V, imax, jmax);
	if(min > restriction){
		min = restriction;
	}

	*dt = tau*min;
}

void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
){
	/* TODO implement */
}


