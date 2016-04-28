#include "uvp.h"
#include <float.h>
#include <math.h>
#include <stdio.h>

/* Funcion to get the discrete Laplacian of any matrix U
at point i,j */
double p_get_laplacian(double u, double u_east, double u_west,
    double u_north, double u_south, double dx_2, double dy_2) {

    return (u_west - 2.0*u + u_east)/dx_2 +
            (u_north - 2.0*u + u_south)/dy_2;
}

/* Calculate the discretized expressions for the momentum equations
according to equations 10 and 11 and apply boundary values
according to equation 18 */
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
    double dx_2 = dx*dx, dy_2 = dy*dy;
    double dU2dx, dUVdy, dUVdx, dV2dy;
    double tmpUx1, tmpUx2, tmpVx1, tmpVx2;
    double tmpUy1, tmpUy2, tmpVy1, tmpVy2;
    double u_ij, v_ij;

    double dx_4 = 1.0/(4.0*dx);
    double dy_4 = 1.0/(4.0*dy);

    /* Calculate F */
    for (j = 1; j <= jmax; j++) {
        for (i = 1; i < imax; i++) {
            u_ij    = U[i][j];

            /* Get discrete Laplacian of U at i,j */
            LapU = p_get_laplacian(u_ij,U[i+1][j],U[i-1][j],
                U[i][j+1],U[i][j-1],dx_2, dy_2);

            /* Compute d(u^2)/dx at i,j */
            tmpUx1  = u_ij + U[i+1][j];
            tmpUx2  = U[i-1][j] + u_ij;
            dU2dx   = ((tmpUx1*tmpUx1 - tmpUx2*tmpUx2) +
                        alpha*(fabs(tmpUx1)*(u_ij - U[i+1][j]) -
                        fabs(tmpUx2)*(U[i-1][j] - u_ij)))*dx_4;

            /* Compute d(u*v)/dx at i,j */
            tmpVx1  = V[i][j] + V[i+1][j];
            tmpVx2  = V[i][j-1] + V[i+1][j-1];
            dUVdy   = ((tmpVx1*(u_ij + U[i][j+1])-tmpVx2*(U[i][j-1] + u_ij)) +
                        alpha*(fabs(tmpVx1)*(u_ij - U[i][j+1]) -
                        fabs(tmpVx2)*(U[i][j-1] - u_ij)))*dy_4;

            /* Compute F at i,j */
            F[i][j] = u_ij + dt*(LapU/Re - dU2dx - dUVdy + GX);
        }
        /* Boundary conditions */
        F[0][j]     = U[0][j];
        F[imax][j]  = U[imax][j];
    }

    /* Calculate G */
    for (i = 1; i <= imax; i++) {
        for (j = 1; j < jmax; j++) {
            v_ij    = V[i][j];

            /* Get discrete Laplacian of V at i,j */
            LapV = p_get_laplacian(v_ij,V[i+1][j],V[i-1][j],
                V[i][j+1],V[i][j-1],dx_2, dy_2);

            /* Compute d(uv)/dy at i,j */
            tmpUy1  = U[i][j] + U[i][j+1];
            tmpUy2  = U[i-1][j] + U[i-1][j+1];
            dUVdx   = ((tmpUy1*(v_ij + V[i+1][j])-tmpUy2*(V[i-1][j] + v_ij)) +
                        alpha*(fabs(tmpUy1)*(v_ij - V[i+1][j]) -
                        fabs(tmpUy2)*(V[i-1][j] - v_ij)))*dx_4;

            /* Compute d(v^2)/dy at i,j */
            tmpVy1  = v_ij + V[i][j+1];
            tmpVy2  = V[i][j-1] + v_ij;
            dV2dy   = ((tmpVy1*tmpVy1 - tmpVy2*tmpVy2) +
                        alpha*(fabs(tmpVy1)*(v_ij - V[i][j+1]) -
                        fabs(tmpVy2)*(V[i][j-1] - v_ij)))*dy_4;

            /* Compute G at i,j */
            G[i][j] = v_ij + dt*(LapV/Re - dUVdx - dV2dy + GY);
        }
        /* Boundary conditions */
        G[i][0]     = V[i][0];
        G[i][jmax]  = V[i][jmax];
    }
}

/* Calculate right hand side of the Pressure Poisson Equation
according to equation 12 */
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
    double dtdx = 1.0/(dt*dx);
    double dtdy = 1.0/(dt*dy);

    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
            RS[i][j] = (F[i][j]-F[i-1][j])*dtdx + (G[i][j]-G[i][j-1])*dtdy;
        }
    }
}

/* Finds the maximum of the absolute values of a matrix */
double p_find_abs_max(double **matrix, int imax, int jmax){
	int i,j;

    /* Initialize maximum value to the first element*/
	double found_max = fabs(matrix[0][0]);

	for(i = 0; i <= imax+1; ++i){
		for(j = 0; j <= jmax+1; ++j){
			if(found_max < fabs(matrix[i][j])){
				found_max = fabs(matrix[i][j]);
			}
		}
	}
	return found_max;
}

/* Calculate the adaptive step size using equation 14 */
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

	restriction = Re/(2.0 * (1.0/(dx*dx) + 1.0/(dy*dy)));
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

/* Calculate the updated U and V
according to equatiosn 8 and 9*/
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
    double dt_dx = dt/dx;
    double dt_dy = dt/dy;

    /* Update U */
	for (i = 1; i < imax; i++) {
	    for (j = 1; j <= jmax; j++) {
	        U[i][j] = F[i][j] - dt_dx*(P[i+1][j]-P[i][j]);
	    }
	}

    /* Update V */
    for (i = 1; i <= imax; i++) {
        for (j = 1; j < jmax; j++) {
            V[i][j] = G[i][j] - dt_dy*(P[i][j+1]-P[i][j]);
        }
    }
}
