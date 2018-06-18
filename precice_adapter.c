#include "precice_adapter.h"
#include "boundary_val.h"
#include <stdlib.h>
#include "adapters/c/SolverInterfaceC.h"
#include <stdio.h>
#include <string.h>

int *precice_set_interface_vertices(char *problem, int imax, int jmax, double dx, double dy, 
									double x_origin, double y_origin,
                                    int num_coupling_cells, int meshID, int **FLAG){
	int dimension   = 3;//precicec_getDimension();
    int* vertexIDs  = (int*)malloc(num_coupling_cells*sizeof(int));

	printf("PROCESS: precice_adapter.c -> Debug_6\n");
	double* vertices = (double*)malloc(num_coupling_cells*dimension*sizeof(double));
    
	/* -------------------CASE 1 GENERALIZED SCAN) BEGINS--------------------- */
	if(strcmp(problem, "F1_heat_exchanger")==0 || strcmp(problem, "F1_heat_exchanger")==0){
		

		int coupledcellcount = 0;
		for(int i = 0; i<imax; i++){
			for(int j=0; j<jmax; j++){ 
				
				// scanning from bottom to top, left to right 
				if(FLAG[i][j]&(1<<9)){

					if(B_N(FLAG[i][j])){
						vertices[dimension*coupledcellcount]     = x_origin + (i - 0.5)*dx;
						vertices[dimension*coupledcellcount + 1] = y_origin + j*dy;
						vertices[dimension*coupledcellcount + 2] = 0;

						coupledcellcount++;
					}

					if(B_S(FLAG[i][j])){
						vertices[dimension*coupledcellcount]     = x_origin + (i - 0.5)*dx;
						vertices[dimension*coupledcellcount + 1] = y_origin + (j - 1)*dy;
						vertices[dimension*coupledcellcount + 2] = 0;

						coupledcellcount++;
					}
				}
			}
		}
		precicec_setMeshVertices(meshID, num_coupling_cells, vertices, vertexIDs);
		/* -------------------CASE 1 GENERALIZED SCAN ENDS--------------------- */
	}

	else{
		/* -------------------------CASE 2 BEGINS------------------------------ */
		printf("PROCESS: precice_adapter.c -> Debug_7\n");
		
		int coupledcellcount = 0;
		for(int j=1; j<jmax-1; j++){ //left boundary
			if(FLAG[0][j]&(1<<9)){
				vertices[dimension*coupledcellcount]     = x_origin;
				vertices[dimension*coupledcellcount + 1] = y_origin + (j - 0.5)*dy;
				vertices[dimension*coupledcellcount + 2] = 0;
				coupledcellcount++;
			}
			
		}
		for(int j=1; j<jmax-1; j++){ //Right boundary 
			if(FLAG[imax-1][j]&(1<<9)){
				vertices[dimension*coupledcellcount]     = x_origin + (imax - 2)*dx;
				vertices[dimension*coupledcellcount + 1] = y_origin + (j - 0.5)*dy;
				vertices[dimension*coupledcellcount + 2] = 0;
				coupledcellcount++;
			}
			
		}
		for(int i=1; i<imax-1; i++){ //Top boundary 
			if(FLAG[i][jmax-1]&(1<<9)){
				vertices[dimension*coupledcellcount]     = x_origin + (i - 0.5)*dx;
				vertices[dimension*coupledcellcount + 1] = y_origin + (jmax - 2)*dy;
				vertices[dimension*coupledcellcount + 2] = 0;
				coupledcellcount++;
			}
			
		}
		for(int i=1; i<imax-1; i++){ //Bottom boundary
			if(FLAG[i][0]&(1<<9)){
				vertices[dimension*coupledcellcount]     = x_origin + (i - 0.5)*dx;
				vertices[dimension*coupledcellcount + 1] = y_origin;
				vertices[dimension*coupledcellcount + 2] = 0;
				coupledcellcount++;
			}
		}
		printf("PROCESS: precice_adapter.c -> Debug_8\n");

		printf("n=%d   c=%d \n", num_coupling_cells,coupledcellcount);
		precicec_setMeshVertices(meshID, num_coupling_cells, vertices, vertexIDs);
	}		
	/* -------------------------CASE 2 ENDS------------------------------ */

    printf("PROCESS: precice_adapter.c -> Debug_9\n");

	free(vertices);
	return vertexIDs;
	
}

void precice_write_temperature(	char *problem, int imax, int jmax, int num_coupling_cells, 
								double *temperature, int *vertexIDs,
                               	int temperatureID, double **TEMP, int **FLAG)
{
	if(strcmp(problem, "F1_heat_exchanger")==0 || strcmp(problem, "F1_heat_exchanger")==0){
	/* ------------------- CASE 1: GENERALIZED CASE: preCICE TEMPERATURE set 
								according to North and South cells --------------------- */	
		int count = 0;
		for(int i = 0; i<imax; i++){
			for(int j=0; j<jmax; j++){ 
				
				// scanning from bottom to top, left to right 
				if(FLAG[i][j]&(1<<9)){

					// Setting preCICE temperature if there is fluid to the north of the coupling cell
					if(B_N(FLAG[i][j])){
						temperature[count] = TEMP[i][j+1];
						count++;
					}
					// Setting preCICE temperature if there is fluid to the south of the coupling cell
					if(B_S(FLAG[i][j])){
						temperature[count] = TEMP[i][j-1];
						count++;
					}
				}
			}
		}
		precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);
	}

	else{
		int count = 0;

		for(int j=1; j<jmax-1; j++){ //left boundary
			if(FLAG[0][j]&(1<<9)){
				temperature[count] = TEMP[1][j];
				count++;
			}
		
		}
		for(int j=1; j<jmax-1; j++){ //Right boundary
			if(FLAG[imax-1][j]&(1<<9)){
				temperature[count] = TEMP[imax-2][j];
				count++;
			}
			
		}
		for(int i=1; i<imax-1; i++){ //Top boundary
			if(FLAG[i][jmax-1]&(1<<9)){
				temperature[count] = TEMP[i][jmax-2];
				count++;
			}
			
		}
		for(int i=1; i<imax-1; i++){ //Bottom boundary
			if(FLAG[i][0]&(1<<9)){
				temperature[count] = TEMP[i][1];
				count++;
			}
			
		}
		precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);
	}
}

void set_coupling_boundary(	char* problem, int imax, int jmax, double dx, double dy, 
							double *heatflux, double **TEMP, int **FLAG)
{
	
/* -------------------CASE 1: GENERALIZED CASE: Domain TEMPERATURE is set according
								to the North and South cells --------------------- */
	if(strcmp(problem, "F1_heat_exchanger")==0 || strcmp(problem, "F1_heat_exchanger")==0){
		int count = 0;
		for(int i = 0; i<imax; i++){
			for(int j=0; j<jmax; j++){ 
				
				
				if(FLAG[i][j]&(1<<9)){
					
					// Setting domain temperature if there is fluid to the north of the coupling cell
					if(B_N(FLAG[i][j])){
					TEMP[i][j]= TEMP[i][j+1]+ dy*(heatflux[count]);
					count++;
					
					}

					// Setting domain temperature if there is fluid to the south of the coupling cell
					if(B_S(FLAG[i][j])){
						TEMP[i][j]= TEMP[i][j-1]+ dy*(heatflux[count]);
						count++;
					}
				}
			}
		}
	}

	else{
		//precice.readBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, heatflux); shift this to main as well
		int count = 0;
		for(int j=1; j<jmax-1; j++){ //left boundary
			if(FLAG[0][j]&(1<<9)){
				TEMP[0][j]= TEMP[1][j]+ dx*(heatflux[count]);
				count++;
			}
			
		}
		for(int j=1; j<jmax-1; j++){ //Right boundary
			if(FLAG[imax-1][j]&(1<<9)){
				TEMP[imax-1][j]= TEMP[imax-2][j]+ dx*(heatflux[count]);
				count++;
			}
			
		}
		for(int i=1; i<imax-1; i++){ //Top boundary
			if(FLAG[i][jmax-1]&(1<<9)){
				TEMP[i][jmax-1]= TEMP[i][jmax-2]+ dy*(heatflux[count]);
				count++;
			}
			
		}
		for(int i=1; i<jmax-1; i++){ //Bottom boundary
			if(FLAG[i][0]&(1<<9)){
				TEMP[i][0]= TEMP[i][1]+ dy*(heatflux[count]);
				count++;
			}
			
		}
	}
}


void write_checkpoint(	double time, double **U, double **V, double **TEMP, 
					double *time_cp, double **U_cp, double **V_cp, double **TEMP_cp, 
					int imax, int jmax){

	time_cp = &time;

	for(int i=0; i<imax; i++){
		for(int j=0; j<jmax; j++){
			U_cp[i][j]		= U[i][j];
			V_cp[i][j]		= V[i][j];
			TEMP_cp[i][j] 	= TEMP[i][j];
		}
	}
}


void restore_checkpoint(double *time, double **U, double **V, double **TEMP,
						double time_cp, double **U_cp, double **V_cp, double **TEMP_cp,
						int imax, int jmax){

	time = &time_cp;

	for(int i=0; i<imax; i++){
		for(int j=0; j<jmax; j++){
			U[i][j]		= U_cp[i][j];
			V[i][j]		= V_cp[i][j];
			TEMP[i][j] 	= TEMP_cp[i][j];
		}
	}	
}

	

	
