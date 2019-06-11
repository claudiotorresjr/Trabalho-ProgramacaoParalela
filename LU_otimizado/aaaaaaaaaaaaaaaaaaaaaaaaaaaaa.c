#include <stdio.h>
#include <math.h>
#include <string.h>

#include <omp.h>

#include "fatLU.h"
#include "immintrin.h"

//#define NTHREADS 1
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

void metodoDeGauss(int threads, double *A, double *b, double *L, int tam)
{
	int j, i, k, blockLin;
	int lin = 0, linend;
	double pivo;
	int NTHREADS = threads;

	memset(L, 0.0, tam*tam*sizeof(double));
	for(i = 0; i < tam; ++i)
	{
		L[i*tam + i] = 1.0;
	}

	for(j = 0; j < tam - 1; ++j)
	{
		//acha todos os m
		pivo = A[j*tam + j];
		#pragma omp parallel for
		for(i = j + 1; i < tam; ++i)
		{
			L[i*tam + j] = A[i*tam + j]/pivo;
			A[i*tam + j] = 0.0;
		}

		#pragma omp parallel default(none) private(blockLin, lin, linend, i, k) \
				shared(NTHREADS, j, A, b, L, tam) num_threads(NTHREADS)
		{
			int ID = omp_get_thread_num(); 
			for(blockLin = ID; blockLin < (tam - j - 1)/BLOCK_SIZE; blockLin += NTHREADS)
			{
				lin = (j + 1) + blockLin*BLOCK_SIZE; linend = lin + BLOCK_SIZE;
				//printf("blockLin %d de %d blocos -- pivo == %d\n", blockLin, (tam - j - 1)/BLOCK_SIZE, j);
				for(i = lin; i < linend; ++i)
				{
					//printf("fazendo elementos de A[%d][%d] ate A[%d][%d]\n\n", i, lin, i, linend-1);
					for(k = j + 1; k < tam; k++)
					{
						A[i*tam + k] -= L[i*tam + j]*A[j*tam + k];
					}
				}
			}
		}
		
		i = (j + 1) + ((tam - j - 1)/BLOCK_SIZE - 1)*BLOCK_SIZE + BLOCK_SIZE;
		if(((tam - j - 1)/BLOCK_SIZE == 0))
		{
			i = j + 1;
		}
		for(; i < tam; ++i)
		{
			for(k = j + 1; k < tam; ++k)
			{
				//printf("Linha %d terminando a coluna %d\n", i, k);
				A[i*tam + k] -= L[i*tam + j]*A[j*tam + k];
			}
			
			//printf("\n");
		}
	}
}