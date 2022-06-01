//CODICE SERIALE

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


// COMMAND LINE:
// cd /mnt/c/Users/black/OneDrive/BackUp/Master\ Archive/Insegnamenti/HPC/Assignment\ B
// gcc densityS.c -o densityS.x -lm
// ./densityS.x 3 2 input.bin
// ./densityS.x 3 2

struct particlesDistribution
{
   	int numParts;
   	float *partsPositions;
};

struct grid
{
	int grid_length;
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
	return ptr;
}

float euclidianDistance(int coord_i, int j_part, struct particlesDistribution pd, int N){
	//from i -> x,y,z
	int coord_z = coord_i % N;
	int coord_y = (coord_i/N) % N;
	int coord_x = (coord_i/(N*N)) % N;
	//compute euclidian distance
	float ed = sqrt((coord_x-(*(pd.partsPositions+j_part))*N)*(coord_x-(*(pd.partsPositions+j_part))*N)+(coord_y-(*(pd.partsPositions+j_part+1))*N)*(coord_y-(*(pd.partsPositions+j_part+1))*N)+(coord_z-(*(pd.partsPositions+j_part+2))*N)*(coord_z-(*((pd.partsPositions+j_part+2))*N)));

	//printf("\ncoord_i %d-> coord_x: %d, coord_y: %d ,coord_z: %d\n", coord_i, coord_x, coord_y, coord_z);
	return ed;
}

float computeDensity(int coord_i, int r, float volume_sfera, struct particlesDistribution pd, int N){
	int nPart=0;
	for(int j=0; j < pd.numParts; j++){ // stampo coordinate delle particelle
		if(euclidianDistance(coord_i, j*3, pd, N) <= r)
			nPart++;
	}
	printf("\nGridPoint-%dnPart:%d\n",coord_i, nPart);
	// calcolare il numero di particelle che cadono nella sfera di raggio R centrata in coord e dividere per il VOLUME della sfera = 4/3*pi*R^3
	return (nPart/volume_sfera);
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
    // INPUT FILE variables
    FILE *fptr;
    struct particlesDistribution pd;
    struct grid g;

    if(argc > 1){

	    // Use current time as seed for random generator
	    srand(time(0));

	    //Create the GRID:
	    g.grid = createGrid(atoi(argv[1]));
	    g.grid_length = atoi(argv[1])*atoi(argv[1])*atoi(argv[1]);

		if(argc == 3) // NO INPUT-FILE!!!
		{
			printf("\nNO INPUT-FILE: let's generate the particlesDistribution randomly...\n");
			//GENERATE THE PARTICLES DISTRIBUTION

			pd.numParts = rand()%g.grid_length; //generated randomly between 0 and N^3
			printf("NumParts:%d\n", pd.numParts);
			pd.partsPositions = (float*)malloc(pd.numParts * 3 * sizeof(float));		

		    for(int i=0; i < pd.numParts*3; i++)
		    	*(pd.partsPositions+i) = (float)rand()/(float)((unsigned)RAND_MAX + 1);
	    }

    	if(argc > 3) // Take the distribution from the INPUT FILE
    	{
    		//READING from INPUT FILE if it doesnt exist than CREATE IT
    		if(fptr = fopen(argv[3],"rb")){
			    // READ the INPUT FILE:
		    	printf("\nREADING...\n");
		    	int nPart;
		    	float coords[3];
		    	// fread(addressData, sizeData, numbersData, pointerToFile);
		    	fread(&nPart, sizeof(nPart), 1, fptr);
		    	printf("Np: %d\n", nPart);

			    //Create the DISTRIBUTION:
			    pd.numParts = nPart;
				pd.partsPositions = (float*)malloc(nPart * 3 * sizeof(float));

		    	for(int i = 0; i < nPart; i++)
				{
			    	fread(&coords, sizeof(coords), 1, fptr);
			    	*(pd.partsPositions+i*3) = coords[0];
			    	*(pd.partsPositions+i*3+1) = coords[1];
			    	*(pd.partsPositions+i*3+2) = coords[2];
			    	printf("x_%d: %f=%f\ty_%d: %f=%f\tz_%d: %f=%f\n", i, coords[0],  *(pd.partsPositions+i*3), i, coords[1], *(pd.partsPositions+i*3+1), i, coords[2], *(pd.partsPositions+i*3+2));
		    	}	

		    	fclose(fptr);	    
  
    		}else{ // input file doesnt exist
    			printf("Input file doesnt exist: let's create it!\n");

			    // CREATE the INPUT FILE:
			    if ((fptr = fopen(argv[3],"wb")) == NULL)
			    {
			        printf("Error! opening file");
			        // Program exits if the file pointer returns NULL.
			        exit(1);
			    }

			    // WRITING in the file:

					//GENERATE RANDOMLY THE PARTICLES DISTRIBUTION
					pd.numParts = rand()%g.grid_length; //generated randomly between 0 and N
				    printf("Np:\t%d\n", pd.numParts); 	

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

		//Compute and Print DENSITIES + Output file:

		//Sphere volume
		float volume_sfera = (4*M_PI*pow(atoi(argv[2]), 3))/3; 
		printf("\nVolume_sfera: %f\n", volume_sfera);

		for(int i=0; i < g.grid_length; i++){
			//computeDensity(coord, raggio, volume_sfera, pd.partsPositions, N)
			*(g.grid+i) = computeDensity(i, atoi(argv[2]), volume_sfera, pd, atoi(argv[1]));
		}
		//print the grid
		printf("\nDensities:");
		for(int i=0; i < g.grid_length; i++){
			printf("\n%d-%f", i, *(g.grid+i));
		}
		printf("\n");

		//DENSITIES OUTPUT-FILE
		fptr = fopen("densities.txt","w");

	    if(fptr == NULL)
	    {
	      	printf("Error!");   
	      	exit(1);             
	   	}

	   	fprintf(fptr,"%d\n\n", atoi(argv[1])); // N 

	   	for(int i=0; i < g.grid_length; i++)
	   	{
	   		fprintf(fptr, "%f\n", *(g.grid+i));
	   	}
	   	fclose(fptr);

		//COMPUTE POTENTIALS + POTENTIALS OUTPUT-FILE
		fptr = fopen("potentials.txt","w");

	    if(fptr == NULL)
	    {
	      	printf("Error!");   
	      	exit(1);             
	   	}

	   	fprintf(fptr,"%d\n\n", pd.numParts); // N 

	   	for(int i=0; i < pd.numParts; i++)
	   	{
	   		//computePotential(particella, distribuzionePart, N)
	   		printf("\nPOTENZIALE TOTALE:\t%f\n", computePotential(i, pd, atoi(argv[1])));
	   		fprintf(fptr, "%f,%f,%f,%f\n", (*pd.partsPositions+i)*atoi(argv[1]), (*pd.partsPositions+i+1)*atoi(argv[1]), (*pd.partsPositions+i+2)*atoi(argv[1]), computePotential(i, pd, atoi(argv[1])));
	   	}
	   	fclose(fptr);
    }
    else{
    	printf("No input from command line!\n");
		return 0;
    }	        
    free(pd.partsPositions);
    free(g.grid);
	return 0;
}