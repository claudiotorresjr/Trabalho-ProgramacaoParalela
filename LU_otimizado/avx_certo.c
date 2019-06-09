#include <stdio.h>
#include <math.h>
#include <string.h>

#include <omp.h>

#include "fatLU.h"
#include "immintrin.h"

#define NTHREADS 8
#define BLOCK_SIZE 4

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

void metodoDeGauss(double *A, double *b, double *L, int tam)
{
	int j, i, k, l, kstart;
	double m;
	double pivo;

	__m256d vetA, vetB, vetC; //vetor que guarda elementos matriz A
	__m256d vetM;
	 //vetor dos coeficientes de multiplicacao (fator de eliminacao)


	memset(L, 0.0, tam*tam*sizeof(double));
	for(i = 0; i < tam; ++i)
	{
		L[i*tam + i] = 1.0;
	}
	
	
	int be = 0;
	for(j = 0; j < tam - 1; ++j)
	{	
		//imprimeMatriz(A, tam);
		//printf("\n");
		//printf("j == %d %d\n", j, be);
		pivo = A[j*tam + j];
		//printf("1 pivo de numero %d e Bloco %d\n", j, be/4);
		vetA = _mm256_loadu_pd(&A[j*tam + be]); //A[j][l]
		//#pragma omp parallel for private(m, vetB, vetM) shared(be, vetA, L, A, b) num_threads(NTHREADS)
		for(i = j + 1; i < tam; ++i)
		{
			//printf("comecando pela linha %d com m == %lf\n", i, A[i*tam + j]/pivo);
			//printf("fazendo emelementos de A[%d][%d] ate A[%d][%d]\n", i, be, i, be+3);
			m = A[i*tam + j]/pivo;
			vetM = _mm256_broadcast_sd(&m);
			vetB = _mm256_loadu_pd(&A[i*tam + be]); //A[i][0]
			vetC = _mm256_mul_pd(vetA, vetM);
			L[i*tam + j] = m;
			vetB = _mm256_sub_pd(vetB, vetC);
			//A[i][l] = ...
			_mm256_storeu_pd(&A[i*tam + be], vetB);
			A[i*tam + j] = 0.0;
			b[i] = b[i] - m*b[j];
		}
		for(k = (j/BLOCK_SIZE)+1; k < tam/BLOCK_SIZE; ++k)
		{	
			kstart = k*BLOCK_SIZE;
			//printf("2 pivo de numero %d e Bloco %d\n", j, k);
			vetA = _mm256_loadu_pd(&A[j*tam + kstart]); //A[j][l]
			{	
				#pragma omp parallel for private(vetB, vetM) shared(kstart, vetA, L, A, b) num_threads(NTHREADS)
				for(i = j + 1; i < tam; ++i)
				{
					//printf("comecando pela linha %d com m == %lf\n", i, A[i*tam + j]/pivo);
					//printf("fazendo emelementos de A[%d][%d] ate A[%d][%d]\n", i, kstart, i, kstart+3);
					vetM = _mm256_broadcast_sd(&L[i*tam + j]);
					vetB = _mm256_loadu_pd(&A[i*tam + kstart]); //A[i][0]
					vetC = _mm256_mul_pd(vetA, vetM);

					vetB = _mm256_sub_pd(vetB, vetC);
					//A[i][l] = ...
					_mm256_storeu_pd(&A[i*tam + kstart], vetB);
				}
			}
		}
		if((j+1)%BLOCK_SIZE == 0 && j != 0)
		{
			be += 4;
		}
	}
}