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

void metodoDeGauss(int threads, double *A, double *b, double *L, int tam)
{
	int j, i, k, blockstart, blockend, blockCol;
	//int lin = 0, linend, blockLin;
	int final;
	double pivo;

	//__m256d vetA, vetB, vetC; //vetor que guarda elementos matriz A
	//__m256d vetM;
	 //vetor dos coeficientes de multiplicacao (fator de eliminacao)
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
		for(i = j + 1; i < tam; ++i)
		{
			L[i*tam + j] = A[i*tam + j]/pivo;
			A[i*tam + j] = 0.0;
		}

		#pragma omp parallel default(none) private(blockCol, blockstart, blockend, i, k) \
				shared(j, A, b, L, tam, NTHREADS)
		{
			int ID = omp_get_thread_num(); 

			//printf("blockLin %d de %d blocos -- pivo == %d\n", blockLin, (tam - j - 1)/BLOCK_SIZE, j);
			for(blockCol = ID; blockCol < (tam - j - 1)/BLOCK_SIZE; blockCol += NTHREADS)
			{	
				blockstart = (j + 1) + blockCol*BLOCK_SIZE; blockend = blockstart + BLOCK_SIZE;
				//printf("blockCol %d de %d blocos -- pivo == %d\n", blockCol, (tam - j - 1)/BLOCK_SIZE, j);
				for(i = j + 1; i < tam; ++i)
				{
					//printf("fazendo elementos de A[%d][%d] ate A[%d][%d]\n\n", i, lin, i, linend-1);
					for(k = blockstart; k < blockend; k++)
					{
						A[i*tam + k] -= L[i*tam + j]*A[j*tam + k];
					}
				}
				//vetB = _mm256_loadu_pd(&A[j*tam + blockstart]);
				//for(i = j + 1; i < tam; ++i)
				//{
				//	//printf("fazendo elementos de A[%d][%d] ate A[%d][%d]\n\n", i, blockstart, i, blockend-1);
				//	vetM = _mm256_broadcast_sd(&L[i*tam + j]);
				//	vetC = _mm256_mul_pd(vetB, vetM);
				//	
				//	vetA = _mm256_loadu_pd(&A[i*tam + blockstart]); 
				//	vetA = _mm256_sub_pd(vetA, vetC);
				//	_mm256_storeu_pd(&A[i*tam + blockstart], vetA);
				//}
			}
		}
		
		for(i = j + 1; i < tam; ++i)
		{
			//final = blockstart + BLOCK_SIZE;
			final = (j + 1) + ((tam - j - 1)/BLOCK_SIZE - 1)*BLOCK_SIZE + BLOCK_SIZE;
			for(k = final; k < tam; ++k)
			{
				//printf("Linha %d terminando a coluna %d\n", i, k);
				A[i*tam + k] -= L[i*tam + j]*A[j*tam + k];
			}
			
			//printf("\n");
			//b[i] = b[i] - m*b[j];
		}
	}
}
