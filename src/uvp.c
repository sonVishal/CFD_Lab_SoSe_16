#include "uvp.h"
#include <float.h>
#include <math.h>
#include <stdio.h>

void get_laplacian(double **U, double dx, double dy,
    int i, int j, double *Lap) {
    *Lap = (U[i+1][j] - 2*U[i][j] + U[i-1][j])/dx/dx +
            (U[i][j+1] - 2*U[i][j] + U[i][j-1])/dy/dy;
}

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
	int i,j;
    double LapU, LapV;
    double dU2dx, dUVdy, dUVdx, dV2dy;
    double tmpUx1, tmpUx2, tmpVx1, tmpVx2;
    double tmpUy1, tmpUy2, tmpVy1, tmpVy2;
    double u_ij, v_ij;

    for (j = 1; j <= jmax; j++) {
        for (i = 1; i < imax; i++) {
            get_laplacian(U, dx, dy, i, j, &LapU);

            u_ij    = U[i][j];
            tmpUx1  = u_ij + U[i+1][j];
            tmpUx2  = U[i-1][j] + u_ij;
            dU2dx   = ((tmpUx1*tmpUx1 - tmpUx2*tmpUx2)+
                        alpha*(fabs(tmpUx1)*(tmpUx1-2*U[i+1][j]) -
                        fabs(tmpUx2)*(tmpUx2-2*u_ij)))/4/dx;

            tmpVx1  = V[i][j] + V[i+1][j];
            tmpVx2  = V[i][j-1] + V[i+1][j-1];
            dUVdy   = ((tmpVx1*(u_ij+U[i][j+1])-tmpVx2*(U[i][j-1]+u_ij)) +
                        alpha*(fabs(tmpVx1)*(u_ij-U[i][j+1]) -
                        fabs(tmpVx2)*(U[i][j-1]-u_ij)))/4/dy;

            F[i][j] = U[i][j] + dt*(LapU/Re - dU2dx - dUVdy + GX);
        }
        /* Boundary conditions */
        F[0][j]     = U[0][j];
        F[imax][j]  = U[imax][j];
    }

    for (i = 1; i <= imax; i++) {
        for (j = 1; j < jmax; j++) {
            get_laplacian(V, dx, dy, i, j, &LapV);

            v_ij    = V[i][j];
            tmpUy1  = U[i][j] + U[i][j+1];
            tmpUy2  = U[i-1][j] + U[i-1][j+1];
            dUVdx   = ((tmpUy1*(v_ij+V[i+1][j])-tmpUy2*(V[i-1][j]+v_ij)) +
                        alpha*(fabs(tmpUy1)*(v_ij-V[i+1][j]) -
                        fabs(tmpUy2)*(V[i-1][j]-v_ij)))/4/dx;

            tmpVy1  = v_ij + V[i][j+1];
            tmpVy2  = V[i][j-1] + v_ij;
            dV2dy   = ((tmpVy1*tmpVy1 - tmpVy2*tmpVy2) +
                        alpha*(fabs(tmpVy1)*(tmpVy1 - 2*V[i][j+1]) -
                        fabs(tmpVy2)*(tmpVy2 - 2*v_ij)))/4/dy;

            G[i][j] = V[i][j] + dt*(LapV/Re - dUVdx - dV2dy + GY);
        }
        /* Boundary conditions */
        G[i][0]     = V[i][0];
        G[i][jmax]  = V[i][jmax];
    }
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
