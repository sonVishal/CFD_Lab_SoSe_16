#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "sor.h"
#include "boundary_val.h"
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

	const char *szFileName;
	const char *szProblem; /* string that describes the output of .vtk files */

	/* Geometry data */
	double  xlength, ylength, dx, dy;
	int     imax, jmax;

	/* Time stepping data */
	double  t = 0, dt, t_end, tau, dt_value, print_tol;
	int 	n = 0;

	/* Pressure iteration data */
	int     it, itermax;
	double	res, eps, omg, alpha;

	/* Problem dependent quantities */
	double  Re, UI, VI, PI, GX, GY;

	/* Solution arrays */
	double	**U, **V, **P, **RS, **F, **G;

	/* Solution dumps for visualization */
	int storageNum = 1;

	if(argn != 2){
		szFileName = "cavity100.dat";
	}else{
		szFileName = args[1];
	}

	szProblem = "pv_files/worksheet1";

	read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength,
			&ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax,
			&eps, &dt_value);

	/* Allocate memory to all the matrices */
	U	= matrix ( 0 , imax+1 , 0 , jmax+1 );
	V	= matrix ( 0 , imax+1 , 0 , jmax+1 );
	P	= matrix ( 0 , imax+1 , 0 , jmax+1 );
	RS	= matrix ( 1 , imax , 1 , jmax );
	F	= matrix ( 0 , imax , 1 , jmax );
	G	= matrix ( 1 , imax , 0 , jmax );

	init_uvp(UI, VI, PI, imax, jmax, U, V, P);

	/* Write the initial conditions for visualization */
	write_vtkFile(szProblem, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
	printf("INFO: write vtk file at time t=%e \t iteration: %i \t file number: %d \n", t, n, storageNum-1);

	/* MAIN LOOP */
	while(t < t_end){
		/* Select dt according to (14) */
		if (tau > 0) {
			calculate_dt(Re, tau,
					&dt, dx, dy, imax, jmax, U, V);
		}

		/* Output of u, v, p values for visualization, if necessary
		   Note the "print_tol" in the condition. This makes it a lot cleaner 
           (see output on console with and without).
           This is done to make the right hand side sligthly greater than the
           time storage value.
		 */
        print_tol = dt*1e-6;
		if(t + dt > storageNum*dt_value + print_tol){
			/* The position of writing the .vtk files in the loop differs from the algorithm
			 * of the assignment sheet. Here we make sure that if 't' is not exactly a multiple
			 * of dt_value a new vtk-file is generated just before 't' gets larger than the
			 * multiple of dt_value.
			 * Example: dt_value = 0.5
			 * t = ... 0.40, 0.49(#), 0.51, 0.55 ... (#) write vtk-file
			 */
			write_vtkFile(szProblem, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
			printf("INFO: write vtk file at time t=%e \t iteration: %i \t file number: %d \n", t, n, storageNum);
			storageNum++;
		}

		/* Set boundary values for u and v according to (15),(16) */
		boundaryvalues(imax,jmax,U,V);

		/* Compute F(n) and G(n) according to (10),(11),(18) */
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy,
				imax, jmax, U, V, F, G);

		/* Compute the right-hand side rs of the pressure equation (12) */
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);

		/* Initialize res to greater than eps so that it enters the sor loop */
		res = 10*eps;

		/* Initialize iteration counter it */
		it = 0;

		/* Perform a SOR iteration according to (19) -- inner loop */
		while (it < itermax && res > eps) {
			sor(omg, dx, dy, imax, jmax, P, RS, &res);
			it++;
		}

        /*Warning when SOR does not converge*/
		if(it >= itermax){
			printf("WARNING: SOR did not converge at time = %f with a residual of %f. "
					"itermax = %i reached. \n", t, res, itermax);
		}

		/* Compute u(n+1) and v(n+1) according to (8),(9) */
		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P);


		t += dt;
		n++;
	}

	/* Output of u, v, p for visualization */
	write_vtkFile(szProblem, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
	printf("INFO: write vtk file at time t=%e \t iteration: %i \t file number: %d \n", t, n, storageNum);

	/*Deallocate matrices*/
	free_matrix (U, 0 , imax+1 , 0 , jmax+1 );
	free_matrix (V, 0 , imax+1 , 0 , jmax+1 );
	free_matrix (P, 0 , imax+1 , 0 , jmax+1 );
	free_matrix (RS, 1 , imax , 1 , jmax );
	free_matrix (F, 0 , imax , 1 , jmax );
	free_matrix (G, 1 , imax , 0 , jmax );

	return 0;
}