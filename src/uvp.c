#include "uvp.h"
#include <float.h>
#include <math.h>


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
	for(i = 0; i<imax; ++i)
		for(j = 0; j < jmax; ++j)
			if(found_max < matrix[i][j])
				found_max = matrix[i][j];

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

	restriction = dx / abs(p_find_max(U, imax, jmax));
	/* TODO: (DL) Think about: maybe just use min function from math.h instead...*/
	if(min > restriction)
		min = restriction;

	restriction = dy / abs(p_find_max(V, imax, jmax));
	if(min > restriction)
		min = restriction;


	/*TODO: (DL) the case when min == 0, is not described in the instructions.
	 * For now it is handled like a negative value (leave at value from read_parameter).
	 */
	if (min > 0){
		*dt = tau*min;
	}else{ /* if(min <= 0) */
		printf("INFO: computed restriction of 'dt' is negative %d. The value remains"
				"at %d", min, *dt);
	}
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

