#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <getopt.h>   /* getopt */

#include <omp.h>

#include "fatLU.h"

static void usage(char *progname)
{
	fprintf(stderr, "Forma de uso: %s [ -n <ordem da matriz> ]\n", progname);
	exit(1);
}

inline double generateRandomA(unsigned int i, unsigned int j, unsigned int k)
{
	double invRandMax = 1.0 / (double)RAND_MAX;
	return ((i==j)?((double)(k<<1)):(1.0))  * ((double)rand() * invRandMax);
}

inline double generateRandomB(unsigned int k)
{
	double invRandMax = 1.0 / (double)RAND_MAX;
	return ((double)(k<<2)) * ((double)rand() * invRandMax);
}

int main(int argc, char *argv[])
{
	int tam;
	int opt;

	srand(20191);

	if(argc != 3)
	{
		usage(argv[0]);
	}
	
	while((opt = getopt(argc, argv, "n:")) != -1)
	{
		switch(opt)
		{
			case 'n':
				tam = atoi(optarg);
				break;
			default:   
				usage(argv[0]);
		}
	}


	
	int i, j;

	double start = omp_get_wtime();
	
	double *b = (double *)malloc(tam*sizeof(double));

	double **A = (double **)malloc(tam*sizeof(double*));		
	for(i = 0; i < tam; i++)
	{
		A[i] = (double *)malloc(tam*sizeof(double));
	}

	double **L = (double **)malloc(tam*sizeof(double*));
	for(i = 0; i < tam; i++)
	{
		L[i] = (double *)malloc(tam*sizeof(double));
	}
	

	for (i = 0; i < tam; ++i) 
	{
		for (j = 0; j < tam; ++j) 
		{	
			A[i][j] = generateRandomA(i, j, tam);
		}
		b[i] = generateRandomB(tam);
	}
	
	double *x = (double *)malloc(tam*sizeof(double));
	double *y = (double *)malloc(tam*sizeof(double));
	
	//LIKWID_MARKER_INIT;
	/*--------------------------------------
	(1)-> A.x = b -> fatorar A em L.U
	(2)-> L.U.x = b
	como U eh triangular, fica facil resolver o sistema:
	(3)-> U.x = y
	tendo y, resolvemos o sistema:
	(4)-> L.y = b
	--------------------------------------*/
	metodoDeGauss(A, b, L, tam);
	
	forwardSubstitution(L, y, b, tam);

	retroSubstitution(A, x, y, tam);
	printf("----------Resultado-----------\n");
	imprimeVetor(x, tam);

	free(A);
	free(L);
	free(b);
	free(x);
	free(y);
	
	double end = omp_get_wtime();
	printf("Time:%f\n", end - start);
	return 0;
}
