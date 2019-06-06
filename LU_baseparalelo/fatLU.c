#include <stdio.h>
#include <math.h>

#include <omp.h>
#include "fatLU.h"

#define NTHREADS 4

void imprimeMatriz(double **A, int tam)
{
	int i, j;
	for(i = 0; i < tam; ++i)
	{
		printf("|");
		for(j = 0; j < tam; ++j)
		{
			printf(" %f", A[i][j]);
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

void forwardSubstitution(double **A, double *x, double *b, int tam)
{	
	int i, j;

	x[0] = b[0]/A[0][0];
	for(i = 1; i < tam; ++i)
	{
		x[i] = b[i];
		for(j = 0; j < i; ++j)
		{
			x[i] -= A[i][j]*x[j];
			
		}
		x[i] /= A[i][i];
	}
}

void retroSubstitution(double **A, double *x, double *b, int tam)
{
	int i, j;

	x[tam - 1] = b[tam - 1]/A[tam - 1][tam - 1];
	for(i = tam - 2; i >= 0; --i)
	{
		x[i] = b[i];
		for(j = tam - 1; j > i; --j)
		{	
			x[i] -= A[i][j]*x[j];
		}
		x[i] /= A[i][i];
	}
}

void trocaLinhas(double **A, double *b, int tam, int k, int l)
{
	int i;
	double aux;

	for(i = 0; i < tam; i++)
	{
		aux = A[k][i];
		A[k][i] = A[l][i];
		A[l][i] = aux;
	}
	aux = b[k];
	b[k] = b[l];
	b[l] = aux;
}

void metodoDeGauss(double **A, double *b, double **L, double **U, int tam)
{
	int j, k, i;

	for(j = 0; j < tam - 1; ++j)
	{
		
		//Pivotamento
		k = j;
		for(i = j + 1; i < tam; ++i)
		{
			if( fabs(A[i][j]) < fabs(A[k][j]))
			{
				k = i;
			}
		}
		trocaLinhas(A, b, tam, k, j);
	}

	for(i = 0; i < tam; ++i)
	{
		L[i][0] = A[i][0];
	}
	for(j = 1; j < tam; ++j)
	{
		U[0][j] = A[0][j]/L[0][0];
	}
	for(i = 0; i < tam; ++i)
	{
		U[i][i] = 1.0;
	}
	
	for(i = 1; i < tam; ++i)
	{
		for(j = 1; j < tam; ++j)
		{
			if(i >= j)
			{
				L[i][j] = A[i][j];
				for(k = 0; k < j; ++k)
				{
					L[i][j] -= L[i][k]*U[k][j];
				}
			}
			else
			{
				U[i][j] = A[i][j];
				for(k = 0; k < j; ++k)
				{
					U[i][j] -= L[i][k]*U[k][j];
				}
				U[i][j] /= L[i][i];
			}
		}
	}
}