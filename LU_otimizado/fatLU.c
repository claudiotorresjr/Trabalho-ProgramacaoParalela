#include <stdio.h>
#include <math.h>
#include <string.h>

#include <omp.h>

#include "fatLU.h"
#include "immintrin.h"

#define NTHREADS 4
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
	int j, i, k, l, lstart;
	double m;
	double pivo;

	__m256d vetA; //vetor que guarda elementos matriz A
	__m256d vetM; //vetor dos coeficientes de multiplicacao (fator de eliminacao)


	memset(L, 0.0, tam*tam*sizeof(double));
	for(i = 0; i < tam; ++i)
	{
		L[i*tam + i] = 1.0;
	}
	
//	pivo = A[0*tam + 0];
//	//printf("pivo de numero %d e Bloco %d\n", j, 0);
//	vetA = _mm256_loadu_pd(&A[0*tam + 0]); //A[j][l]
//	for(i = 0 + 1; i < tam; ++i)
//	{
//		//printf("comecando pela linha %d\n", i);
//		m = A[i*tam + 0]/pivo;
//		vetM = _mm256_broadcast_sd(&m);
//		vetA = _mm256_mul_pd(vetA, vetM);
//		L[i*tam + 0] = m;
//		for(l = 0 + 1; l < 4; ++l)
//		{
//			//printf("atualizando valores de A[%d][%d] a A[%d][%d]\n", i, l, i, l+3);
//			A[i*tam + l] -= vetA[l];
//		}
//		A[i*tam + 0] = 0.0;
//		b[i] = b[i] - m*b[0];
//	}
//	lstart = l;

	for(j = 0; j < tam - 1; ++j)
	{	
		pivo = A[j*tam + j];

		//if(j < tam - NTHREADS + 1)
		//{

			for(k = (j/BLOCK_SIZE); k < tam/BLOCK_SIZE; ++k)
			{
				//printf("pivo de numero %d e Bloco %d\n", j, k);
				vetA = _mm256_loadu_pd(&A[j*tam + lstart]); //A[j][l]
				//#pragma omp parallel default(none) private(m, i, vetM, vetA, vetA, l) \
				shared(A, L, b, tam, j, pivo) num_threads(NTHREADS) 
				{	
					//#pragma omp for
					for(i = j + 1; i < tam; ++i)
					{
						//printf("comecando pela linha %d\n", i);
						m = A[i*tam + j]/pivo;
						vetM = _mm256_broadcast_sd(&m);
						vetA = _mm256_mul_pd(vetA, vetM);
						L[i*tam + j] = m;
						for(l = j + 1; l < 4; ++l)
						{
							//printf("atualizando valores de A[%d][%d] a A[%d][%d]\n", i, l, i, l+3);
							A[i*tam + l] -= vetA[l];
						}
						A[i*tam + j] = 0.0;
						b[i] = b[i] - m*b[j];
					}
					lstart = l;
				}
			}
		//}
		//else
		//{
		//	for(i = j + 1; i < tam; ++i)
		//	{
		//		m = A[i*tam + j]/pivo;
		//		
		//		L[i*tam + j] = m;
		//		vetM = _mm256_broadcast_sd(&m);
		//		A[i*tam + j] = 0.0;
		//		for(l = j + 1; tam - l >= 4; l+=4)
		//		{
		//			//carrego os 256bits nao alinhados começando em (A[i*tam + j]) para o vetor AVX vetA
		//			//faco o mesmo com (A[j*tam + l])
		//			//a cada iteraçao pego os proximos 256bits 
		//			vetA = _mm256_loadu_pd(&A[i*tam + l]); //A[i][l]
		//			vetA = _mm256_loadu_pd(&A[j*tam + l]); //A[j][l]
		//			//A[i][l] - m*A[j][l];
		//			vetA = _mm256_mul_pd(vetA, vetM);
		//			vetA = _mm256_sub_pd(vetA, vetA);
		//			//A[i][l] = ...
		//			_mm256_storeu_pd(&A[i*tam + l], vetA);
		//		}
		//		for (; l < tam; ++l)
		//		{
		//			A[i*tam + l] = A[i*tam + l] - m*A[j*tam + l];
		//		}
		//		b[i] = b[i] - m*b[j];
		//	}
		//}
	}
}