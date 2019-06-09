#include <stdio.h>
#include <omp.h>

#define NTHREADS 64

#include "matriz.h"

void escreveMatriz(double *C, unsigned int tam)
{
	unsigned int i, j;
	for(i = 0; i < tam; ++i)
	{
		for(j = 0; j < tam; ++j)
		{
			printf("%lf ", C[i*tam + j]);
		}
		printf("\n");
	}
}

void multMatrizNormal(double *A, double *B, double *C, unsigned int tam)
{
	unsigned int i, j, k;

	for(i = 0; i < tam; ++i)
	{
		for(j = 0; j < tam; ++j)
		{
			for(k = 0; k < tam; ++k)
			{
				C[i*tam + j] += A[i*tam + k]*B[k*tam + j];
			}
		}
	}
	//escreveMatriz(C, tam);	
}

void multMatrizTransposta(double *A, double *Bt, double *C, unsigned int tam)
{
	unsigned int i, j, k;

	#pragma omp parallel for private(j, k) num_threads(NTHREADS)
	for(i = 0; i < tam; ++i)
	{
		for(j = 0; j < tam; ++j)
		{
			for(k = 0; k < tam; ++k)
			{
				C[i*tam + j] += A[i*tam + k]*Bt[j*tam + k];
			}
		}
	}
	//escreveMatriz(C, tam);	
}

void multMatrizNormalBloco(double *A, double *B, double *C, unsigned int tam)
{
	unsigned int m, r, p;
	unsigned int i, istart, iend;
	unsigned int j, jstart, jend;
	unsigned int k, kstart, kend;

	for (i = 0; i < tam/BLOCK_SIZE; ++i)
		{
			istart = i*BLOCK_SIZE; iend = istart+BLOCK_SIZE; 
			for (j = 0; j < tam/BLOCK_SIZE; ++j) 
			{
				jstart = j*BLOCK_SIZE; jend = jstart+BLOCK_SIZE;
				for (k = 0; k < tam/BLOCK_SIZE; ++k)
				{	
					kstart = k*BLOCK_SIZE; kend = kstart+BLOCK_SIZE;
					for(m = istart; m < iend; ++m)
					{
						for(r = jstart; r < jend; r += Unroll)
						{	
							for(p = kstart; p < kend; ++p)
							{
								C[m*tam + r]   += A[m*tam + p]*B[p*tam + r];
								C[m*tam + r+1] += A[m*tam + p]*B[p*tam + r+1];
								C[m*tam + r+2] += A[m*tam + p]*B[p*tam + r+2];
								C[m*tam + r+3] += A[m*tam + p]*B[p*tam + r+3];
								C[m*tam + r+4] += A[m*tam + p]*B[p*tam + r+4];
								C[m*tam + r+5] += A[m*tam + p]*B[p*tam + r+5];
								C[m*tam + r+6] += A[m*tam + p]*B[p*tam + r+6];
								C[m*tam + r+7] += A[m*tam + p]*B[p*tam + r+7];
							}
						}
					}
				}
			}
		}
	//escreveMatriz(C, tam);	
}
void multMatrizTranspostaBloco(double *A, double *Bt, double *C, unsigned int tam)
{
	unsigned int m, r, p;
	unsigned int i, istart, iend;
	unsigned int j, jstart, jend;
	unsigned int k, kstart, kend;

	for (i = 0; i < tam/BLOCK_SIZE; ++i)
		{
			istart = i*BLOCK_SIZE; iend = istart+BLOCK_SIZE; 
			for (j = 0; j < tam/BLOCK_SIZE; ++j) 
			{
				jstart = j*BLOCK_SIZE; jend = jstart+BLOCK_SIZE;
				for (k = 0; k < tam/BLOCK_SIZE; ++k)
				{	
					kstart = k*BLOCK_SIZE; kend = kstart+BLOCK_SIZE;
					for(m = istart; m < iend; ++m)
					{
						for(r = jstart; r < jend; r += Unroll)
						{	
							for(p = kstart; p < kend; ++p)
							{
								C[m*tam + r]   += A[m*tam + p]*Bt[r*tam + p];
								C[m*tam + r+1] += A[m*tam + p]*Bt[(r+1)*tam + p];
								C[m*tam + r+2] += A[m*tam + p]*Bt[(r+2)*tam + p];
								C[m*tam + r+3] += A[m*tam + p]*Bt[(r+3)*tam + p];
								C[m*tam + r+4] += A[m*tam + p]*Bt[(r+4)*tam + p];
								C[m*tam + r+5] += A[m*tam + p]*Bt[(r+5)*tam + p];
								C[m*tam + r+6] += A[m*tam + p]*Bt[(r+6)*tam + p];
								C[m*tam + r+7] += A[m*tam + p]*Bt[(r+7)*tam + p];
							}
						}
					}
				}
			}
		}
	//escreveMatriz(C, tam);	
}