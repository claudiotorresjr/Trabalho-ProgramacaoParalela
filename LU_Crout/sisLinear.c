#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <getopt.h>   /* getopt */

#include <omp.h>

#include "fatLU.h"

static void usage(char *progname)
{
	fprintf(stderr, "Forma de uso: %s [ -p <numero de threads> -n <ordem da matriz> ]\n", progname);
	exit(1);
}

double generateRandomA(unsigned int i, unsigned int j, unsigned int k)
{
	double invRandMax = 1.0 / (double)RAND_MAX;
	return ((i==j)?((double)(k<<1)):(1.0))  * ((double)rand() * invRandMax);
}

double generateRandomB(unsigned int k)
{
	double invRandMax = 1.0 / (double)RAND_MAX;
	return ((double)(k<<2)) * ((double)rand() * invRandMax);
}

int main(int argc, char *argv[])
{
	int tam, threads;
	int opt;

	srand(20191);

	if(argc != 5)
	{
		usage(argv[0]);
	}
	
	while((opt = getopt(argc, argv, "p:n:")) != -1)
	{
		switch(opt)
		{
			case 'n':
				tam = atoi(optarg);
				break;
			case 'p':
				threads = atoi(optarg);
				break;
			default:   
				usage(argv[0]);
		}
	}

	double start = omp_get_wtime();

	omp_set_num_threads(threads);
	
	int i, j;
	double *A = (double*)aligned_alloc(64, tam*tam*sizeof(double));	
	double *L = (double*)aligned_alloc(64, tam*tam*sizeof(double));
	double *U = (double*)aligned_alloc(64, tam*tam*sizeof(double));	
	
	double *b = (double*)aligned_alloc(64, tam*sizeof(double));
	double *x = (double*)aligned_alloc(64, tam*sizeof(double));
	double *y = (double*)aligned_alloc(64, tam*sizeof(double));

	for (i = 0; i < tam; ++i) 
	{
		for (j = 0; j < tam; ++j) 
		{	
			A[i*tam + j] = generateRandomA(i, j, tam);
		}
		b[i] = generateRandomB(tam);
	}

	/*--------------------------------------
	(1)-> A.x = b -> fatorar A em L.U
	(2)-> L.U.x = b
	como U eh triangular, fica facil resolver o sistema:
	(3)-> U.x = y
	tendo y, resolvemos o sistema:
	(4)-> L.y = b
	--------------------------------------*/

	metodoDeGauss(A, b, L, U, tam);
	
	forwardSubstitution(L, y, b, tam);

	retroSubstitution(U, x, y, tam);
	printf("----------Resultado-----------\n");
	imprimeVetor(x, tam);

	free(A);
	free(L);
	free(U);
	free(b);
	free(x);
	free(y);

	double end = omp_get_wtime();
	printf("Time:%f\n", end - start);

	return 0;
}