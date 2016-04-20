#include "uvp.h"
#include <float.h>
#include <math.h>
#include <stdio.h>

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
    int i,j;

    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
            RS[i][j] = ((F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy)/dt;
        }
    }
}

double p_find_abs_max(double **matrix, int imax, int jmax){
	int i,j;

	/* equals: minimum (negative) double value; note the "-" (!!)
	 * https://stackoverflow.com/questions/1153548/minimum-double-value-in-c-c*/
	double foundMax = -DBL_MAX;

	/*TODO: (DL) this is most likely wrong, check out again how the matrix is structured */
	for(i = 0; i<imax; ++i){
		for(j = 0; j < jmax; ++j){
			if(foundMax < fabs(matrix[i][j])){
				foundMax = fabs(matrix[i][j]);
			}
		}
	}
	return foundMax;
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

	restriction = dx / p_find_abs_max(U, imax, jmax);
	min = fmin(min, restriction);

	restriction = dy / p_find_abs_max(V, imax, jmax);
	min = fmin(min, restriction);

	*dt = tau*min;

	if(*dt <= 0){
		printf("WARNING: dt has invalid value: %e", *dt);
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
    int i,j;

	for (i = 1; i < imax; i++) {
	    for (j = 1; j < jmax; j++) {
	        U[i][j] = F[i][j] - dt*(P[i+1][j]-P[i][j])/dx;
            V[i][j] = G[i][j] - dt*(P[i][j+1]-P[i][j])/dy;
	    }
        U[i][jmax] = F[i][jmax] - dt*(P[i+1][jmax]-P[i][jmax])/dx;
	}

    for (j = 1; j < jmax; j++) {
        V[imax][j] = G[imax][j] - dt*(P[imax][j+1]-P[imax][j])/dy;
    }
}
