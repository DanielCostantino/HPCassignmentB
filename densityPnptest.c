// CODICE PARALLELO

#define _GNU_SOURCE
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#include <sched.h>
#include <omp.h>
#include "mpi.h"

// gcc version 7.4.0 (Ubuntu 7.4.0-1ubuntu1~18.04)

// ssh ct1-005.area.trieste.it -l costa
// github_token: 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
// github_token2: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

// COMMAND LINE:
// echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
// cd /mnt/c/Users/black/Documents/GitHub/HPCassignmentB
// mpicc -fopenmp densityPnptest.c -o densityPnptest.o -lm
// mpirun -np 1 densityPnptest.o
// mpirun -np 1 densityPnptest.o 2 2 input.bin
// mpirun -np 8 densityPnptest.o 2 2 input.bin
// mpirun -np 8 densityPnptest.o 2 2

// mpirun --mca pml ob1 --mca btl tcp,self -np 8 densityPnptest.o 2 2 input.bin  //pml -> ob1; btl -> tcp
// mpirun --mca pml ucx --mca btl tcp,self --map-by node -np 8 densityPnptest.o 2 2 //pml -> ucx; btl -> tcp

int rank;

struct particlesDistribution
{
   	int numParts;
   	float *partsPositions;
};

struct grid
{
	int numGrid_points;
	float *grid;
};

float* createGrid(int N){
	//create the 3d-grid
	float* ptr = (float*)malloc(N * N * N * sizeof(float));
    if (ptr == NULL)
    {
        fprintf(stderr, "Out of memory!");
        exit(0);
    }

	/*/printGrid
	printf("pow(N,3)= %f\n", pow(N,3));
	for (int i = 0; i < pow(N,3); ++i)
		printf("%f,", *(ptr+i));
	printf("\n");*/

	return ptr;
}

float euclidianDistance(int coord_i, int j_part, struct particlesDistribution pd, int N, int gridPoint){
	//from i -> x,y,z
	int coord_z = rank + coord_i % gridPoint;
	int coord_y = rank + (coord_i/gridPoint) % gridPoint;
	int coord_x = rank + (coord_i/(gridPoint*gridPoint)) % gridPoint;
	//compute euclidian distance
	float ed = sqrt((coord_x-(*(pd.partsPositions+j_part))*N)*(coord_x-(*(pd.partsPositions+j_part))*N)+(coord_y-(*(pd.partsPositions+j_part+1))*N)*(coord_y-(*(pd.partsPositions+j_part+1))*N)+(coord_z-(*(pd.partsPositions+j_part+2))*N)*(coord_z-(*((pd.partsPositions+j_part+2))*N)));

	//printf("\ncoord_i %d-> coord_x: %d, coord_y: %d ,coord_z: %d\n", coord_i, coord_x, coord_y, coord_z);
	return ed;
}

float computeDensity(int coord_i, int r, float volume_sfera, struct particlesDistribution pd, int N, int gridPoint){
	int nPart=0;
	for(int j=0; j < pd.numParts; j++){ // stampo coordinate delle particelle
		if(euclidianDistance(coord_i, j*3, pd, N, gridPoint) <= r)
			nPart++;
	}
	//printf("\nGridPoint-%d has nPart:%d\n",coord_i, nPart);
	// calcolare il numero di particelle che cadono nella sfera di raggio R centrata in coord e dividere per il VOLUME della sfera = 4/3*pi*R^3
	return nPart/volume_sfera;
}

float eD(int i, int j, struct particlesDistribution pd, int N){
	return sqrt((((*(pd.partsPositions+i))*N)-(*(pd.partsPositions+j))*N)*(((*(pd.partsPositions+i))*N)-(*(pd.partsPositions+j))*N)+(((*(pd.partsPositions+i+1))*N)-(*(pd.partsPositions+j+1))*N)*(((*(pd.partsPositions+i+1))*N)-(*(pd.partsPositions+j+1))*N)+(((*(pd.partsPositions+i+2))*N)-(*(pd.partsPositions+j+2))*N)*(((*(pd.partsPositions+i+2))*N)-(*((pd.partsPositions+j+2))*N)));
}

float computePotential(int i_part, struct particlesDistribution pd, int N){
	float p = 0;
	float ed = 0;

	for(int i=0; i < pd.numParts; i++)
		if(i!=i_part)
			p+= 1/eD(i, i_part, pd, N);
		//printf("\nDistanza:\t%f",eD(i, i_part, pd, N));		
		//printf("\nPart-%d con potenziale p_%d:\t%f", i_part, i, ed);				
	return p;
}

int main(int argc, char *argv[])
{
	// PERFORMANCE variables
	double start_time, end_time, t;  

    // INPUT FILE variables
    FILE *fptr;

    struct particlesDistribution pd;
    struct grid g;
    float* grid_output = NULL;
    
	// In MPI tutte le variabili sono private del processo (vengono creati dei linux processes, che hanno sempre memoria privata!!)
    int	numprocs;
    int ROOT = 0;

    int len;
    char hostname[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv); // initializeMPI

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs); // tells to MPIprocesses "HOW MUCH we are"
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // tells to each MPIprocess its own rank

    MPI_Status status;

	int block_dim;
	

    if(argc > 1){

   	MPI_Get_processor_name(hostname, &len);
//	   printf ("Number of tasks= %d My rank= %d Running on %s\n", numprocs,rank,hostname);

		if(rank==ROOT)
		{	
    		start_time = MPI_Wtime();
    	}
	    // Use current time as seed for random generator
	    srand(time(0));
		// split grid in "np blocks" and assign each of them to one MPIprocess
		//TO DO: round-robin assignment if block_dim is not integer
		block_dim = pow(atoi(argv[1]), 3) / numprocs;
		//printf("\nblock_dim:\t%d\n", block_dim);

		//each MPIprocess(root included) create an own block-grid of dimension block_dim:
	    g.grid = createGrid(block_dim);
	    g.numGrid_points = block_dim;

	    //ROOT generates/reads the "particles distribution"
		if(rank == ROOT)
		{		
			if(argc == 3) // NO INPUT-FILE!!!
			{
				//printf("\nNO INPUT-FILE: let's generate the particlesDistribution randomly...\n");
				//GENERATE THE PARTICLES DISTRIBUTION
				pd.numParts = atoi(argv[2]); // <--- NpTest | rand()%(int)(pow(atoi(argv[1]), 3)); //generated randomly between 0 and N^3
				//printf("NumParts:%d\n", pd.numParts);
				pd.partsPositions = (float*)malloc(pd.numParts * 3 * sizeof(float));  // ?? sto rompendo questo malloc qua ??		

			    for(int i=0; i < pd.numParts*3; i++){					
			    	*(pd.partsPositions+i) = ((float)rand())/((float)((unsigned)RAND_MAX + 1));
			    	//printf("*(pd.partsPositions+i)=%f\n", *(pd.partsPositions+i));
				}
			    /*/print particles distribution
			    for(int i=0; i < pd.numParts; i++)
					printf("Part-%d: %f, %f, %f\n", i, *(pd.partsPositions+i*3), *(pd.partsPositions+i*3+1), *(pd.partsPositions+i*3+2));
				*/
		    }

	    	if(argc > 3) // Take the distribution from the INPUT FILE
	    	{
	    		//READING from INPUT FILE if it doesnt exist than CREATE IT
	    		if(fptr = fopen(argv[3],"rb")){
				    // READ the INPUT FILE:
			    	//printf("\nREADING...\n");
			    	int nPart;
			    	float coords[3];
			    	// fread(addressData, sizeData, numbersData, pointerToFile);
			    	fread(&nPart, sizeof(nPart), 1, fptr);
			    	//printf("Np: %d\n", nPart);

				    //Create the DISTRIBUTION:
				    pd.numParts = nPart;
					pd.partsPositions = (float*)malloc(nPart * 3 * sizeof(float));

			    	for(int i = 0; i < nPart; i++)
					{
				    	fread(&coords, sizeof(coords), 1, fptr);
				    	*(pd.partsPositions+i*3) = coords[0];
				    	*(pd.partsPositions+i*3+1) = coords[1];
				    	*(pd.partsPositions+i*3+2) = coords[2];
				    	//printf("x_%d: %f=%f\ty_%d: %f=%f\tz_%d: %f=%f\n", i, coords[0],  *(pd.partsPositions+i*3), i, coords[1], *(pd.partsPositions+i*3+1), i, coords[2], *(pd.partsPositions+i*3+2));
			    	}	

			    	fclose(fptr);	    
	  
	    		}else{ // input file doesnt exist
	    			printf("Input file doesnt exist: let's create it!\n");

				    // CREATE the INPUT FILE:
				    if ((fptr = fopen(argv[3],"wb")) == NULL)
				    {
				        //printf("Error! opening file");
				        // Program exits if the file pointer returns NULL.
				        exit(1);
				    }
				    // WRITING in the file:
						//GENERATE RANDOMLY THE PARTICLES DISTRIBUTION
						pd.numParts = rand()%g.numGrid_points; //generated randomly between 0 and N
					    //printf("Np:\t%d\n", pd.numParts); 	

					    //Create the DISTRIBUTION:
						pd.partsPositions = (float*)malloc(pd.numParts * 3 * sizeof(float));

					    for(int i=0; i< pd.numParts*3;i++)
					    	*(pd.partsPositions+i) = (float)rand()/(float)((unsigned)RAND_MAX + 1);

				    fwrite(&pd.numParts, sizeof(pd.numParts), 1, fptr); 
				    for(int i = 0; i < pd.numParts*3; i++)
						fwrite(pd.partsPositions+i, sizeof(float), 1, fptr);

				    fclose(fptr); 	   
	    		}		
	    	}
	    }
    		
	    // ONPENMPI & DENSITIES:
		    //sending to each MPIprocess the particles distribution
			MPI_Bcast(&pd.numParts, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
			//printf("Rank: %d, NumParts:%d\n", rank, pd.numParts);
			//printf("pd.partsPositions=%p\n",pd.partsPositions);
		    if(rank!=ROOT) // each MPIprocess != ROOT 
		    {
		    	//printf("rank-%dpd.numParts=%d",rank,pd.numParts);
				pd.partsPositions = (float*)malloc(pd.numParts*3*sizeof(float));
				if(pd.partsPositions==NULL){
					printf("Rank:%d-pd.partsPositions==NULL\n", rank);		
				    exit(1); 
				}		
		    }
			MPI_Bcast(pd.partsPositions, pd.numParts*3, MPI_FLOAT, ROOT, MPI_COMM_WORLD);
			//for(int i=0; i < pd.numParts; i++)
			//	printf("Part-%d: %f, %f, %f\n", i, *(pd.partsPositions+i*3), *(pd.partsPositions+i*3+1), *(pd.partsPositions+i*3+2));

			//compute DENSITIES
			//Sphere volume

			float radius = atoi(argv[1])*sqrt(2*pow(atoi(argv[1]),2))*(1.0/cbrt(atoi(argv[2]))); 
			float volume_sfera = (4*M_PI*pow(radius,3))/3; 
			//printf("\nVolume_sfera: %f\n", volume_sfera);
			if(rank==ROOT){
				printf("radius: %f\n", radius);
				printf("volume_sfera: %f\n", volume_sfera);
			}
			//printf("g.numGrid_points=%d\n", g.numGrid_points);
			for(int i=0; i < g.numGrid_points; i++)
				*(g.grid+i) = computeDensity(i, radius, volume_sfera, pd, atoi(argv[1]), g.numGrid_points);
							//computeDensity(coord, raggio, volume_sfera, pd.partsPositions, N, g.numGrid_points)

			/*/print the grid
			printf("\nDensities:");
			for(int i=0; i < g.numGrid_points; i++)
				printf("\n%d-%f", i, *(g.grid+i));

			printf("\n");*/			
				  			     	
			//DENSITIES OUTPUT-FILE
			if(rank == ROOT){
				grid_output = createGrid(atoi(argv[1]));
				
				if (grid_output == NULL)
				{
			      	printf("grid_output, ciaone error!");   
			      	exit(1);             
			   	}				
			}
		
			MPI_Gather(g.grid, g.numGrid_points, MPI_FLOAT, grid_output, g.numGrid_points, MPI_FLOAT, ROOT, MPI_COMM_WORLD);
						
			if(rank == ROOT)
			{
				fptr = fopen("densities.txt","w");

			    if(fptr == NULL)
			    {
			      	printf("Error!");   
			      	exit(1);             
			   	}
		
			   	fprintf(fptr,"%d\n\n", atoi(argv[1])); // N 

			   	for(int i=0; i < pow(atoi(argv[1]),3); i++)
			   		fprintf(fptr, "%f\n", *(grid_output+i));
			   	
			   	fclose(fptr);
			 }
		
			// OPENMP + POTENTIALS
			if(rank == ROOT){
				#pragma omp parallel num_threads(8) proc_bind(close)
				{
			    /*#pragma omp master
			    {
			      char *proc_bind_names[] = { "false (no binding)",
							  "true",
							  "master",
							  "close",
							  "spread" };
			      
			      
			      // get the current binding
			      int binding = omp_get_proc_bind();

				//printf(" proc bind is set to \"%s\"\n", proc_bind_names[binding] );
			    }

				int my_thread_id = omp_get_thread_num();  // private variable	                                           									      			    		
				int cpu_num = sched_getcpu();
			    int place   = omp_get_place_num();
			
				//printf("Place is %d, my CPU is %d, my omp thread is %d, my mpi thread is %d\n" ,place, cpu_num, my_thread_id, rank);
		    	*/
				#pragma omp single		
					fptr = fopen("potentials.txt","w");
			
				    	if(fptr == NULL)
				    	{
				      		printf("Error!");   
				      		exit(1);             
				   		}

						#pragma omp single
				   			fprintf(fptr,"%d\n\n", pd.numParts); // Np
						
				   		#pragma omp for schedule(static)
				   		for(int i = 0; i < pd.numParts; i++)
				   		{
				   			//computePotential(particella, distribuzionePart, N)
				   			//printf("\nPOTENZIALE TOTALE:\t%f\n", computePotential(i, pd, atoi(argv[1])));
				   			//printf("interation_%d taken from thread-%d\n", i, my_thread_id);
				   			fprintf(fptr, "%f,%f,%f,%f\n", (*pd.partsPositions+i)*atoi(argv[1]), (*pd.partsPositions+i+1)*atoi(argv[1]), (*pd.partsPositions+i+2)*atoi(argv[1]), computePotential(i, pd, atoi(argv[1])));
				   		}

				   		#pragma omp single
				   			fclose(fptr);
	
    				//printf( "\tGreetings from thread num %d\n", my_thread_id);			
				}
			}

			//timing performance
			if(rank==ROOT)
			{			
				end_time = MPI_Wtime();
				t = end_time - start_time;
				printf ( "%10.8f ", t);
			}

	}
    else{
    	printf("No input from command line!\n");
	    MPI_Finalize();	
		return 0;
    }

    // FREE the dynamic arrays, remember that u cannot free without a malloc/calloc!	
    MPI_Barrier(MPI_COMM_WORLD);
    if(pd.partsPositions!=NULL)	
    	free(pd.partsPositions);	
    if(g.grid!=NULL)
    	free(g.grid);    	
    if(grid_output!=NULL)
    	free(grid_output); 
    	   
    MPI_Finalize();
	return 0;
}

