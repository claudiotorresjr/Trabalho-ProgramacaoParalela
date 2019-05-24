#include <stdio.h>
#include <math.h>

#include "fatLU.h"


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
	printf("%d\n", tam-1);
	for(i = tam - 2; i >= 0; --i)
	{
		//printf("%d %d\n", i, j);
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

void metodoDeGauss(double **A, double *b, double **L, int tam)
{
	int j, k, i, l;
	double m;

	for(i = 0; i < tam; i++)
	{
		L[i][i] = 1.0;
	}
	
	for(j = 0; j < tam - 1; ++j)
	{
		
		//Pivotamento
		k = j;
		for(i = j + 1; i < tam; ++i)
		{
			if( fabs(A[i][j]) > fabs(A[k][j]))
			{
				k = i;
			}
		}
		trocaLinhas(A, b, tam, k, j);
		
		for(i = j + 1; i < tam; ++i)
		{		
			m = A[i][j]/A[j][j];
			L[i][j] = m;
			A[i][j] = 0.0;

			for(l = j+1; l < tam; l++)
			{
				A[i][l] = A[i][l] - m*A[j][l];
			}
			b[i] = b[i] - m*b[j];
		}
	}
}

/*void fatoracaoLU(double **A, double **L, int tam)
{
	int j, i, l;
	double m;
	
	for(i = 0; i < N; i++)
	{
		L[i][i] = 1.0;
	}

	for(j = 0; j < n - 1; ++j)
	{
		for(i = j + 1; i < n; ++i){		
			m = A[i][j]/A[j][j];
			L[i][j] = m;
			A[i][j] = 0;

			for(l = j + 1; l < n; ++l)
			{
				A[i][l] = A[i][l] - m*A[j][l];
			}
		}		
	}
} */