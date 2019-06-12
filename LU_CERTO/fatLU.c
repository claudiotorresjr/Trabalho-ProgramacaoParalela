#include <stdio.h>
#include <math.h>
#include <string.h>

#include <omp.h>

#include "fatLU.h"
#include "immintrin.h"

//#define NTHREADS 8
#define BLOCK_SIZE 8

void imprimeMatriz(double *A, int tam)
{
	int i, j;
	for(i = 0; i < tam; ++i)
	{
		printf("|");
		for(j = 0; j < tam; ++j)
		{
			printf(" %f", A[i*tam + j]);
		}
		printf(" |\n");
	}
}

void imprimeVetor(double *v, int tam)
{
	int i;
	for(i = 0; i < tam; ++i)
	{
		printf("|%f|\n",v[i]);
	}
}

void forwardSubstitution(double *A, double *x, double *b, int tam)
{	
	int i, j;

	x[0] = b[0]/A[0];
	for(i = 1; i < tam; ++i)
	{
		x[i] = b[i];
		for(j = 0; j < i; ++j)
		{
			x[i] -= A[i*tam + j]*x[j];
			
		}
		x[i] /= A[i*tam + i];
	}
}

void retroSubstitution(double *A, double *x, double *b, int tam)
{
	int i, j;

	x[tam - 1] = b[tam - 1]/A[(tam - 1)*tam + (tam - 1)];
	for(i = tam - 2; i >= 0; --i)
	{
		x[i] = b[i];
		for(j = tam - 1; j > i; --j)
		{	
			x[i] -= A[i*tam + j]*x[j];
		}
		x[i] /= A[i*tam + i];
	}
}

void trocaLinhas(double *A, double *b, int tam, int k, int l)
{
	int i;
	double aux;

	for(i = 0; i < tam; i++)
	{
		aux = A[k*tam + i];
		A[k*tam + i] = A[l*tam + i];
		A[l*tam + i] = aux;
	}
	aux = b[k];
	b[k] = b[l];
	b[l] = aux;
}

void metodoDeGauss(double *A, double *b, double *L, double *U, int tam)
{
	int j, i, k, m;
	int blockstart, blockend, blockCol;
	//int lin = 0, linend, blockLin;
	int final;
	double pivo;
	double *Ut = (double*)aligned_alloc(64, tam*tam*sizeof(double));	
	
	__m256d vetA, vetB, vetC; //vetor que guarda elementos matriz A
	__m256d vetM;
	 //vetor dos coeficientes de multiplicacao (fator de eliminacao)
	int NTHREADS = omp_get_num_threads();

	memset(L, 0.0, tam*tam*sizeof(double));
	memset(Ut, 0.0, tam*tam*sizeof(double));

	unsigned int cont = 0;
	for(i = 0; i < tam; ++i)
	{
		L[i*tam + 0] = A[i*tam + 0];
		Ut[i*tam + i] = 1;
	}
	for(j = 1; j < tam; ++j)
	{
		Ut[j*tam + 0] = A[0*tam + j]/L[0*tam + 0];
	}

	for(m = 1; m < tam; ++m)
	{
		
		for(i = m; i < tam; ++i)
		{
		
			for(k = 0; k < m; ++k)
			{
		
				A[i*tam + m] -= L[i*tam + k]*Ut[m*tam + k];
			}
			L[i*tam + m] = A[i*tam + m];
		}

		for(j = m + 1; j < tam;  ++j)
		{
		
			for(k = 0; k < m;  ++k)
			{	
				A[m*tam + j] -= L[m*tam + k]*Ut[j*tam + k];
			}
			Ut[j*tam + k] = A[m*tam + j]/L[m*tam + m];
		}
	}


	for(i = 0; i < tam; ++i)
	{
		for(j = 0; j < tam; ++j)
		{
			U[i*tam + j] = Ut[j*tam + i];
		}
	}
}
