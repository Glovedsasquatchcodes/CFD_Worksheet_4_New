#include "helper.h"
#include "visual.h"
#include "init.h"
#include"uvp.h"
#include"boundary_val.h"
#include"sor.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "string.h"
#include "adapters/c/SolverInterfaceC.h"
#include "precice_adapter.h"
#include "adapters/c/Constants.h"


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
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
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
int main(int argc, char* argv[]){
			char *problem = argv[1];
			char *filename = malloc(strlen(problem) + 5);
			strcpy(filename, problem);
	                strcat(filename, ".dat");

			//define parameter variables
			double Re;                /* reynolds number   */
			double UI;                /* velocity x-direction */
			double VI;                /* velocity y-direction */
			double PI;                /* pressure */
			double GX;                /* gravitation x-direction */
			double GY;                /* gravitation y-direction */
			double t_end;             /* end time */
			double xlength;           /* length of the domain x-dir.*/
			double ylength;           /* length of the domain y-dir.*/
			double dt;                /* time step */
			double dx;                /* length of a cell x-dir. */
			double dy;                /* length of a cell y-dir. */
			int  imax;                /* number of cells x-direction*/
			int  jmax;                /* number of cells y-direction*/
			double alpha;             /* uppwind differencing factor*/
			double omg;               /* relaxation factor */
			double tau;               /* safety factor for time step*/
			int  itermax;             /* max. number of iterations for pressure per time step  */
			double eps;               /* accuracy bound for pressure*/
			double dt_value;           /* time for output */
			double Pr;
			double TI;
			//double T_h;
			//double T_c;
			double beta;
			double x_origin;
		        double y_origin;
		        char *geometry;
		        char *precice_config;
		        char *participant_name;
		        char *mesh_name;
		        char *read_data_name;
		        char *write_data_name;

			//Read and assign the parameter values from file
			read_parameters(filename, &imax, &jmax, &xlength, &ylength, &dt, &t_end, &tau, &dt_value, &eps, &omg, &alpha, &itermax,&GX, &GY, &Re, &Pr, &UI, &VI, &PI, &TI, &beta, &dx, &dy, &x_origin, &y_origin, geometry, precice_config, participant_name, mesh_name, read_data_name, write_data_name);

			geometry = "natural_convection.pgm";			
			printf("Name of geometry: %s\n", geometry);
			precice_config = "precice-configs/precice_config_plate_implicit.xml";
			participant_name = "Fluid";
			mesh_name = "Fluid-Mesh";
			read_data_name = "Heat-Flux";
			write_data_name = "Temperature";

			//Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap
			printf("Allocate the matrices for P(pressure), U(velocity_x), V(velocity_y), F, and G on heap... \n");
			double **P = matrix(0, imax-1, 0, jmax-1);
			double **U = matrix(0, imax-1, 0, jmax-1);
			double **V = matrix(0, imax-1, 0, jmax-1);
			double **F = matrix(0, imax-1, 0, jmax-1);
			double **G = matrix(0, imax-1, 0, jmax-1);
			double **RS = matrix(0, imax-1, 0, jmax-1);
			int **flag = imatrix(0, imax-1, 0, jmax-1);
			double **U_cp = matrix(0, imax-1, 0, jmax-1);
			double **V_cp = matrix(0, imax-1, 0, jmax-1);

			double **T;
			T = matrix(0, imax-1, 0, jmax-1);


			double **T_cp;
			T_cp = matrix(0, imax-1, 0, jmax-1);

			int num_coupling;	
			printf("Matrices allocated... \n \n");

			//Flag Initialization
			init_flag(geometry, imax, jmax, flag, &num_coupling);
			
			//Initialize the U, V and P
			init_uvp(UI, VI, PI, TI, imax, jmax, U, V, P, T, flag);


			//Make solution folder
			struct stat st = {0};
			char sol_folder[80];
			sprintf( sol_folder,"Results_%s",problem);
			if (stat(sol_folder, &st) == -1) {
		    		mkdir(sol_folder, 0700);
			}

			char sol_directory[80];
			sprintf( sol_directory,"Results_%s/sol", problem);
			
			// Algorithm starts from here
			printf("Alogrithm started........\n");
			double t=0;
			double time_cp=0;
			int n=0;
			int n1=0;
			printf("Debug_1\n");

			// initialize preCICE
			precicec_createSolverInterface(participant_name, precice_config, 0, 1);
			printf("Debug_2\n");

			const char* coric = precicec_actionReadIterationCheckpoint(); 
			const char* cowic = precicec_actionWriteIterationCheckpoint();

			int dim = precicec_getDimensions();

			printf("Debug_3\n");
			
			// define coupling mesh
			int meshID = precicec_getMeshID(mesh_name);
			printf("Debug_4\n");
			//int num_coupling_cells = num_coupling;//determine no. of coupling cells 
			printf("Debug_5\n");					
    			int* vertexIDs = precice_set_interface_vertices(imax,jmax, dx, dy, x_origin, y_origin, num_coupling, meshID, flag); // get coupling cell ids
			printf("Debug_11\n");

			// define Dirichlet part of coupling written by this solver
			int temperatureID = precicec_getDataID(write_data_name, meshID);
			double* temperatureCoupled = (double*) malloc(sizeof(double) * num_coupling);

			printf("Debug_12\n");

			// define Neumann part of coupling read by this solver
			int heatFluxID = precicec_getDataID(read_data_name, meshID);
			double* heatfluxCoupled = (double*) malloc(sizeof(double) * num_coupling);
				
			// call precicec_initialize()
			double precice_dt = precicec_initialize();

			// initialize data at coupling interface
			precice_write_temperature(imax,jmax,num_coupling,temperatureCoupled,vertexIDs,temperatureID,T,flag); 
			precicec_initialize_data(); // synchronize with OpenFOAM
			precicec_readBlockScalarData(heatFluxID, num_coupling, vertexIDs, heatfluxCoupled);//read heatfluxCoupled	changed vertexsize to num_coupling_cells	



			while (precicec_isCouplingOngoing()) {
				
				if(precicec_isActionRequired(cowic)){
					printf("\n \n \n Setting the Check point \n \n \n \n");
					write_checkpoint(t, U, V, T, &time_cp, U_cp, V_cp, T_cp, imax, jmax); // save checkpoint
					precicec_fulfilledAction(cowic);
 				}
	
				calculate_dt(Re,tau,&dt,dx,dy,imax,jmax, U, V, Pr);
				dt = fmin(dt, precice_dt); // change solver_dt to dt

				spec_boundary_val(imax, jmax, U, V, flag);

				boundaryvalues(imax, jmax, U, V, flag);

				set_coupling_boundary(imax, jmax, dx, dy, heatfluxCoupled, T, flag); 							
				calculate_temp(T, Pr, Re, imax, jmax, dx, dy, dt, alpha, U, V, flag, TI);
				
				calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G,flag, beta, T);
									
				calculate_rs(dt,dx,dy,imax,jmax,F,G,RS,flag);
											
				int it = 0;
				double res = 10.0;

				do{
					sor(omg,dx,dy,imax,jmax,P,RS,&res,flag);
						++it;
				} while(it<itermax && res>eps);


				printf("SOR iterations = %d ,residual = %f \n", it-1, res);

				if((it==itermax)&&(res>eps)){
					printf("WARNING: SOR Iteration limit reached before convergence. \n");
				}

				calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P,flag);
			
				precice_write_temperature(imax,jmax,num_coupling,temperatureCoupled,vertexIDs,temperatureID,T, flag);
				precice_dt = precicec_advance(dt); // advance coupling
				precicec_readBlockScalarData(heatFluxID, num_coupling, vertexIDs, heatfluxCoupled); //changed vertexsize to num_coupling_cells	
				if(precicec_isActionRequired(coric)){ // timestep not converged
					printf("\n \n \n Restoring the Check point\n \n \n ");
					restore_checkpoint(&t, U, V, T, time_cp, U_cp, V_cp, T_cp, imax, jmax); // set variables back to checkpoint
					precicec_fulfilledAction(coric);
				}
				else
				{
				t =t+ dt;
				n = n+ 1;
				}

				reset_obstacles(U, V, P, T, flag, imax, jmax);
	
				if ((t >= n1*dt_value)&&(t!=0.0)){
					write_vtkFile(sol_directory ,n ,xlength ,ylength ,imax-2 ,jmax-2,dx ,dy ,U ,V ,P,T, x_origin, y_origin);
					printf("Writing Solutions at %f seconds in the file \n",n1*dt_value);
				    	n1=n1+ 1;
				    	continue;
				}

			}


			printf("Algorithm successfully executed...\n \n");
			printf("Freeing dynamically allocated memory...\n");
			    //Free memory
			free_matrix( P, 0, imax-1, 0, jmax-1);
			free_matrix( U, 0, imax-1, 0, jmax-1);
			free_matrix( V, 0, imax-1, 0, jmax-1);
			free_matrix( F, 0, imax-1, 0, jmax-1);
			free_matrix( G, 0, imax-1, 0, jmax-1);
			free_matrix(RS, 0, imax-1, 0, jmax-1);
			free_imatrix(flag, 0, imax-1, 0, jmax-1);
			free_matrix(T, 0, imax-1, 0, jmax-1);
			
			

			printf("End \n");
			return -1;
    
}
