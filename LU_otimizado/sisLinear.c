#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <getopt.h>   /* getopt */

#include <omp.h>
#include <likwid.h>

#include "fatLU.h"

static void usage(char *progname)
{
	fprintf(stderr, "Forma de uso: %s [ -n <ordem da matriz> ]\n", progname);
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
	double *A = (double*)malloc(tam*tam*sizeof(double));	
	double *L = (double*)malloc(tam*tam*sizeof(double));

	double *b = (double*)malloc(tam*sizeof(double));
	double *x = (double*)malloc(tam*sizeof(double));
	double *y = (double*)malloc(tam*sizeof(double));

	for (i = 0; i < tam; ++i) 
	{
		for (j = 0; j < tam; ++j) 
		{	
			A[i*tam + j] = generateRandomA(i, j, tam);
		}
		b[i] = generateRandomB(tam);
	}
	
	//LIKWID_MARKER_INIT;
	/*--------------------------------------
	(1)-> A.x = b -> fatorar A em L.U
	(2)-> L.U.x = b
	como U eh triangular, fica facil resolver o sistema:
	(3)-> U.x = y
	tendo y, resolvemos o sistema:
	(4)-> L.y = b
	--------------------------------------*/

	//puts("----------Matriz A-----------");
	//imprimeMatriz(A, tam);
	//puts("\n");
	//puts("----------Vetor de coeficientes b-----------");
	//imprimeVetor(b, tam);
	//puts("\n");
	//imprimeMatriz(A, tam);
	//puts("----------LU-----------");
	//fatoracaoLU(A,L,tam);
	//LIKWID_MARKER_START("fatLU");
	metodoDeGauss(A, b, L, tam);
	//LIKWID_MARKER_STOP("fatLU");
	//puts("----------Vetor b apos Gauss-----------");
	//imprimeVetor(b, tam);
	//puts("-------------------------");
	//puts("Apos Gauss:");
	//puts("U:");
	//imprimeMatriz(A, tam);
	//puts("\nL:");
	//imprimeMatriz(L, tam);
	
	forwardSubstitution(L, y, b, tam);
	//puts("y:");
	//imprimeVetor(y, tam);
	
	//apos Gauss, A virou U
	retroSubstitution(A, x, y, tam);
	puts("----------Resultado-----------");
	imprimeVetor(x, tam);
	/*
	imprimeMatriz(A);
	imprimeVetor(b);
	
	printf("\n Resultado: \n");
	imprimeVetor(x);
	*/

	free(A);
	free(L);
	free(b);
	free(x);
	free(y);

	//LIKWID_MARKER_CLOSE;
	return 0;
}