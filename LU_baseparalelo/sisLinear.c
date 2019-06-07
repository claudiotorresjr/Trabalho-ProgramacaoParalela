#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <getopt.h>   /* getopt */

//#include <omp.h>
//#include <likwid.h>

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
		memset(L[i], 0, tam);
	}

	double **U = (double **)malloc(tam*sizeof(double*));
	for(i = 0; i < tam; i++)
	{
		U[i] = (double *)malloc(tam*sizeof(double));
		memset(U[i], 0, tam);
	}
	

	//for (i = 0; i < tam; ++i) 
	//{
	//	for (j = 0; j < tam; ++j) 
	//	{	
	//		A[i][j] = generateRandomA(i, j, tam);
	//	}
	//	b[i] = generateRandomB(tam);
	//}

	A[0][0] = 9.0; A[0][1] = 3.0; A[0][2] = 3.0; A[0][3] = 3.0;
	A[1][0] = 3.0; A[1][1] = 10.0; A[1][2] = -2.0; A[1][3] = -2.0;
	A[2][0] = 3.0; A[2][1] = -2.0; A[2][2] = 18.0; A[2][3] = 10.0;
	A[3][0] = 3.0; A[3][1] = -2.0; A[3][2] = 10.0; A[3][3] = 10.0;

	//U[0][0] = 1.0; U[0][1] = 0.333333; U[0][2] = 0.333333; U[0][3] = 0.333333;
	//U[1][0] = 0.0; U[1][1] = 1.0; U[1][2] = -0.333333; U[1][3] = -0.333333;
	//U[2][0] = 0.0; U[2][1] = 0.0; U[2][2] = 1.0; U[2][3] = 1/2;
	//U[3][0] = 0.0; U[3][1] = 0.0; U[3][2] = 0.0; U[3][3] = 1.0;

	// L[0][0] = 9.0; L[0][1] = 0.0; L[0][2] = 0.0; L[0][3] = 0.0;
	// L[1][0] = 3.0; L[1][1] = 9.0; L[1][2] = 0.0; L[1][3] = 0.0;
	// L[2][0] = 3.0; L[2][1] = -3.0; L[2][2] = 16.0; L[2][3] = 0.0;
	// L[3][0] = 3.0; L[3][1] = -3.0; L[3][2] = 8.0; L[3][3] = 4.0;


	b[0] = 24.0; b[1] = 17.0; b[2] = 45.0; b[3] = 49.0;
	
	double *x = (double *)malloc(tam*sizeof(double));
	double *y = (double *)malloc(tam*sizeof(double));
	
	puts("----------Matriz A-----------");
	imprimeMatriz(A, tam);
	puts("\n");	

	metodoDeGauss(A, b, L, U, tam);
	
	puts("U:");
	imprimeMatriz(U, tam);
	
	puts("\nL:");
	imprimeMatriz(L, tam);

	double **R = (double **)malloc(tam*sizeof(double*));
	for(i = 0; i < tam; i++)
	{
		R[i] = (double *)malloc(tam*sizeof(double));
		memset(R[i], 0, tam);
	}
	
	for (i = 0; i < tam; ++i) {
		for (j = 0; j < tam; ++j) {
			for (int k = 0; k < tam; ++k) {
				R[i][j] += L[i][k]*U[k][j];
			}
		}
	}

	puts("\nR:");
	imprimeMatriz(R, tam);

	forwardSubstitution(L, y, b, tam);
	retroSubstitution(U, x, y, tam);

	puts("----------Resultado-----------");
	imprimeVetor(x, tam);

	for (i = 0 ; i < tam ; ++i) {
		free(A[i]);
		free(L[i]);
		free(U[i]);
		free(R[i]);
	}
	free(A);
	free(R);
	free(L);
	free(U);
	free(b);
	free(x);
	free(y);
	
	//LIKWID_MARKER_CLOSE;
	return 0;
}
