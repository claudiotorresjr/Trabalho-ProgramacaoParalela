#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <getopt.h>   /* getopt */

#include <likwid.h>

#include "matriz.h"

static void usage(char *progname)
{
	fprintf(stderr, "Forma de uso: %s [ -n <ordem da matriz> ]\n", progname);
	exit(1);
}

int main(int argc, char *argv[])
{
	int tam, opt;

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

	double *A, *B, *Bt, *C;
	unsigned i, j;

	A =  (double *)aligned_alloc(64, tam*tam*sizeof(double));
	B =  (double *)aligned_alloc(64, tam*tam*sizeof(double));
	Bt = (double *)aligned_alloc(64, tam*tam*sizeof(double));
	C =  (double *)aligned_alloc(64, tam*tam*sizeof(double));

	for (i = 0; i < tam; ++i) 
	{
		for (j = 0; j < tam; ++j) 
		{	
			A[i*tam + j] = i/(j+1.0);
			B[i*tam + j] = i-j;
			C[i*tam + j] = 0.0;
		}
	}

	for (i = 0; i < tam; ++i) 
	{
		for (j = 0; j < tam; ++j) 
		{	
			Bt[i*tam + j] = B[j*tam + i];
		}
	}

	LIKWID_MARKER_INIT;

	LIKWID_MARKER_START("multMatrizNormal");
	multMatrizNormal(A, B, C, tam);
	LIKWID_MARKER_STOP("multMatrizNormal");

	LIKWID_MARKER_START("multMatrizTransposta");
	multMatrizTransposta(A, Bt, C, tam);
	LIKWID_MARKER_STOP("multMatrizTransposta");

	LIKWID_MARKER_START("multMatrizNormalBloco");
	multMatrizNormalBloco(A, B, C, tam);
	LIKWID_MARKER_STOP("multMatrizNormalBloco");

	LIKWID_MARKER_START("multMatrizTranspostaBloco");
	multMatrizTranspostaBloco(A, Bt, C, tam);
	LIKWID_MARKER_STOP("multMatrizTranspostaBloco");

	LIKWID_MARKER_CLOSE;

	return 0;	
}