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

int main(int argc, char *argv[])
{
	int tam;
	int opt;

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
	}
	

	//for (i = 0; i < tam; ++i) 
	//{
	//	for (j = 0; j < tam; ++j) 
	//	{	
	//		A[i][j] = i/(j+1.0);
	//	}
	//}
	//for (i = 0; i < tam; ++i) 
	//{	
	//	b[i] = i*3/(j+4.0);
	//}

	A[0][0] = 3.0; A[0][1] = 2.0; A[0][2] =  4.0;  b[0] = 4.0;
	A[1][0] = 1.0; A[1][1] = 1.0; A[1][2] =  2.0;  b[1] = 2.0;
	A[2][0] = 4.0; A[2][1] = 3.0; A[2][2] = -2.0;  b[2] = 3.0;
	
	double *x = (double *)malloc(tam*sizeof(double));
	double *y = (double *)malloc(tam*sizeof(double));
	
	/*--------------------------------------
	(1)-> A.x = b -> fatorar A em L.U
	(2)-> L.U.x = b
	como U eh triangular, fica facil resolver o sistema:
	(3)-> U.x = y
	tendo y, resolvemos o sistema:
	(4)-> L.y = b
	--------------------------------------*/

	puts("----------LU-----------");
	//fatoracaoLU(A,L,tam);
	metodoDeGauss(A, b, L, tam);
	puts("-------------------------");
	puts("Apos Gauss:");
	puts("U:");
	imprimeMatriz(A, tam);
	puts("\nL:");
	imprimeMatriz(L, tam);
	
	forwardSubstitution(L, y, b, tam);
	//puts("y:");
	//imprimeVetor(y, tam);
	
	//apos Gauss, A virou U
	retroSubstitution(A, x, y, tam);
	puts("Resultado:");
	imprimeVetor(x, tam);
	/*
	imprimeMatriz(A);
	imprimeVetor(b);
	
	printf("\n Resultado: \n");
	imprimeVetor(x);
	*/
	return 0;
}