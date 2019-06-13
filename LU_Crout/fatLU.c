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

	double *Ut = (double*)aligned_alloc(64, tam*tam*sizeof(double));
	
	__m256d vetL, vetU, result; //vetor que guarda elementos matriz A
	__m256d vetS;

	memset(L, 0.0, tam*tam*sizeof(double));
	memset(Ut, 0.0, tam*tam*sizeof(double));
	
	//int NTHREADS = omp_get_num_threads();

	for(i = 0; i < tam; ++i)
	{
		L[i*tam + 0] = A[i*tam + 0];
		Ut[i*tam + i] = 1;
	}
	for(j = 1; j < tam; ++j)
	{
		Ut[j*tam + 0] = A[0*tam + j]/L[0*tam + 0];
	}


	for(m = 1; m < BLOCK_SIZE; ++m)
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

			for(k = 0; k < m; ++k)
			{	
				A[m*tam + j] -= L[m*tam + k]*Ut[j*tam + k];
			}
			Ut[j*tam + m] = A[m*tam + j]/L[m*tam + m];
		}
	}
	for(; m < tam; ++m)
	{	
		#pragma omp parallel default(none) \
			private(i, k, j, vetL, vetU, vetS, result) \
			shared(m, tam, A, L, Ut)
		{
			//int ID = omp_get_thread_num();
			#pragma omp for
			for(i = m; i < tam; ++i)
			{	

				vetS = _mm256_setzero_pd();
				for(k = 0; k + 4 < m; k += 4)
				{	

					vetL = _mm256_loadu_pd(&L[i*tam + k]);
					vetU = _mm256_loadu_pd(&Ut[m*tam + k]);
					//multiplica os valores do vetL com o vetU, e soma com cada um de vetS
					//resultado é salvo em vetS
					vetS = _mm256_fmadd_pd(vetL, vetU, vetS);
					//A[i*tam + m] = A[i*tam + m] - L[i*tam + k]*Ut[m*tam + k];
				}
				//soma os elementos de cada vetor dois a dois alternadamente. cada um possui 4 elementos:
				//a1-a2-a3-a4 e b1-b2-b3-b4 -> result = (a1+a2)-(b1+b2)-(a3+a4)-(b3+b4)
				result = _mm256_hadd_pd(vetS, vetS);
				//como queremos o somatorio dos valores de vetS e nao duas vezes ele,
				//pegamos o result[0] + result[2] --> (a1+a2) + (a3+a4) == somatorio de vetS (a1+a2+a3+a4)
				//mudamos result de __m256d para doble
				A[i*tam + m] -= ((double*)&result)[0] + ((double*)&result)[2];
				for(; k < m; ++k)
				{

					A[i*tam + m] -= L[i*tam + k]*Ut[m*tam + k];
				}
				L[i*tam + m] = A[i*tam + m];
			}
			//#pragma omp barrier

			#pragma omp for
			for(j = m + 1; j < tam; ++j)
			{
				vetS = _mm256_setzero_pd();
				for(k = 0; k + 4 < m; k += 4)
				{	

					vetL = _mm256_loadu_pd(&L[m*tam + k]);
					vetU = _mm256_loadu_pd(&Ut[j*tam + k]);
					vetS = _mm256_fmadd_pd(vetL, vetU, vetS);
				}
				result = _mm256_hadd_pd(vetS, vetS);
				
				A[m*tam + j] -= ((double*)&result)[0] + ((double*)&result)[2];
				for(; k < m; ++k)
				{

					A[m*tam + j] -= L[m*tam + k]*Ut[j*tam + k];
				}
				Ut[j*tam + m] = A[m*tam + j]/L[m*tam + m];
			}
		}
	}
	

	//transposta de Ut (U)
	for(i = 0; i < tam; ++i)
	{
		for(j = 0; j < tam; ++j)
		{
			U[i*tam + j] = Ut[j*tam + i];
		}
	}

	free(Ut);
}

//https://nptel.ac.in/courses/111107062/module2/lecture3/lecture3.pdf
