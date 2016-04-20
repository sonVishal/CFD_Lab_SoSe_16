#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "sor.h"
#include <stdio.h>


    /**
    * The main operation reads the configuration file, initializes the scenario and
    * contains the main loop. So here are the individual steps of the algorithm:
    *
    * - read the program configuration file using read_parameters()
    * - set up the matrices (arrays) needed using the matrix() command
    * - create the initial setup init_uvp(), init_flag(), output_uvp()
    * - perform the main loop
    * - trailer: destroy memory allocated and do some statistics
    *
    * The layout of the grid is described by the first figure below, the enumeration
    * of the whole grid is given by the second figure. All the unknowns correspond
    * to a two dimensional degree of freedom layout, so they are not stored in
    * arrays, but in a matrix.
    *
    * @image html grid.jpg
    *
    * @image html whole-grid.jpg
    *
    * Within the main loop the following big steps are done (for some of the
    * operations a definition is defined already within uvp.h):
    *
    * - calculate_dt() Determine the maximal time step size.
    * - boundaryvalues() Set the boundary values for the next time step.
    * - calculate_fg() Determine the values of F and G (diffusion and confection).
    *   This is the right hand side of the pressure equation and used later on for
    *   the time step transition.
    * - calculate_rs()
    * - Iterate the pressure poisson equation until the residual becomes smaller
    *   than eps or the maximal number of iterations is performed. Within the
    *   iteration loop the operation sor() is used.
    * - calculate_uv() Calculate the velocity at the next time step.
    */
int main(int argn, char** args){

 /**************************
  * Initializing variables *
  *************************/

  /* TODO: Consider making this less annoying by placing function parameters
   *       on a line so one does not have to scroll a long time to look at
   *       something. Right now I am using automatic folds, but I am not sure
   *       if Benjamin has anything to automize this.
   */

	const char *szFileName;
    double  Re, UI;
    double  VI;
    double  PI;
    double  GX;
    double  GY;
    double  t = 0;
    double  t_end;
    double  xlength;
    double  ylength;
    double  dt;
    double  dx;
    double  dy;
    int     imax;
    int     jmax;
    double  alpha;
    double  omg;
    double  tau;
    int     itermax;
    double  eps;
    double  dt_value;
    int 	n = 0;

    double** U;
    double** V;
    double** P;

	if(argn != 2){
		szFileName = "cavity100.dat";
	}else{
		szFileName = args[1];
	}

    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
    V = matrix ( 0 , imax+1 , 0 , jmax+1 );
    P = matrix ( 0 , imax+1 , 0 , jmax+1 );

read_parameters(
    szFileName,
    &Re,
    &UI,
    &VI,
    &PI,
    &GX,
    &GY,
    &t_end,
    &xlength,
    &ylength,
    &dt,
    &dx,
    &dy,
    &imax,
    &jmax,
    &alpha,
    &omg,
    &tau,
    &itermax,
    &eps,
    &dt_value
);

init_uvp(
    UI,
    VI,
    PI,
    imax,
    jmax,
    U,
    V,
    P
);

/* MAIN LOOP */
while(t < t_end){
	/* Select δt according to (14) */
	calculate_dt( Re, tau,
	    &dt, dx, dy, imax, jmax, U, V );

	/* Set boundary values for u and v according to (15),(16) */

	/* Compute F(n) and G(n) according to (10),(11),(18) */

	/* Compute the right-hand side rs of the pressure equation (12) */

	/* Perform a SOR iteration according to (19) -- inner loop */

	/* Compute u(n+1) and v(n+1) according to (8),(9) */

	/* Output of u, v, p values for visualization, if necessary */

	t += dt;
	n++;
}

/* Output of u, v, p for visualization */


/*
 * TODO reminder: deallocation of matrices etc.
 */

/*Deallocate matrices*/ 
    free_matrix (U, 0 , imax+1 , 0 , jmax+1 );
    free_matrix (V, 0 , imax+1 , 0 , jmax+1 );
    free_matrix (P, 0 , imax+1 , 0 , jmax+1 );
return 0;
}
