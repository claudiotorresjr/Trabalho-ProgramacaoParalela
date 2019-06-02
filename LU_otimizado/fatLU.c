#include <stdio.h>
#include <math.h>
#include <string.h>

#include <omp.h>

#include "fatLU.h"
#include "immintrin.h"

#define NTHREADS 2

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
	int j, i, k, l;
	double pivo;

	__m256d vetA; //vetor que guarda elementos matriz A
	__m256d vetB; //vetor nao sei pq criei
	__m256d vetM; //vetor dos coeficientes de multiplicacao (fator de eliminacao)


	memset(L, 0.0, tam*tam*sizeof(double));
	for(i = 0; i < tam; ++i)
	{
		L[i*tam + i] = 1.0;
	}
	
	for(j = 0; j < tam - 1; ++j)
	{	
		//Pivotamento
		//k = j;
		//for(i = j + 1; i < tam; ++i)
		//{
		//	if( ABS(A[i*tam + j]) > ABS(A[k*tam + j]))
		//	{
		//		k = i;
		//	}
		//}
		//trocaLinhas(A, b, tam, k, j);

		pivo = A[j*tam + j];

		//if(j < tam - NTHREADS + 1)
		//{
			#pragma omp parallel default(none) private(i, vetM, vetA, vetB, l) \
			shared(A, L, b, tam, j) num_threads(NTHREADS) 
			{	

				int ID = omp_get_thread_num() + 1;
				double m;
				//#pragma omp for private(m, vetM, vetA, vetB, l)
				for(i = j + ID; i < tam; i+=NTHREADS)
				{			
					m = A[i*tam + j]/A[j*tam + j];
					printf("th %d pegando a linha %d de pivo %d (o m sera dado pela pos: (%d, %d)%lf / (%d, %d)%lf == %lf\n", \
												ID-1, i, j, i,j, A[i*tam + j], j,j, A[j*tam + j], m);
					//printf("th %d calculou m == %lf\n", ID-1, m);
					L[i*tam + j] = m;
					vetM = _mm256_broadcast_sd(&m);

					A[i*tam + j] = 0.0;
					for(l = j + 1; l < tam - (tam % 4); l+=4)
					{
						//carrego os 256bits nao alinhados começando em (A[i*tam + j]) para o vetor AVX vetA
						//faco o mesmo com (A[j*tam + l])
						//a cada iteraçao pego os proximos 256bits 
						vetA = _mm256_loadu_pd(&A[i*tam + l]); //A[i][l]
						vetB = _mm256_loadu_pd(&A[j*tam + l]); //A[j][l]
						//A[i][l] - m*A[j][l];
						vetB = _mm256_mul_pd(vetB, vetM);
						vetA = _mm256_sub_pd(vetA, vetB);
						//A[i][l] = ...
						_mm256_storeu_pd(&A[i*tam + l], vetA);
					}
					for (; l < tam; ++l)
					{
						A[i*tam + l] = A[i*tam + l] - m*A[j*tam + l];
					}
					b[i] = b[i] - m*b[j];
				}
			}
			imprimeMatriz(A, tam);
		//}
		//else
		//{
		//	for(i = j + 1; i < tam; ++i)
		//	{			
		//		//printf("th %d finalizando a linha %d de pivo %d\n", 0, i, j);
		//		m = A[i*tam + j]/pivo;
		//		L[i*tam + j] = m;
		//		vetM = _mm256_broadcast_sd(&m);
		//		
		//		A[i*tam + j] = 0.0;
		//		for(l = j + 1; l < tam - (tam % 4); l+=4)
		//		{
		//			//carrego os 256bits nao alinhados começando em (A[i*tam + j]) para o vetor AVX vetA
		//			//faco o mesmo com (A[j*tam + l])
		//			//a cada iteraçao pego os proximos 256bits 
		//			vetA = _mm256_loadu_pd(&A[i*tam + l]); //A[i][l]
		//			vetB = _mm256_loadu_pd(&A[j*tam + l]); //A[j][l]
		//			//A[i][l] - m*A[j][l];
		//			vetB = _mm256_mul_pd(vetB, vetM);
		//			vetA = _mm256_sub_pd(vetA, vetB);
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