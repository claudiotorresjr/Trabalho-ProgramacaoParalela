#ifndef _FATLU_H
#define _FATLU_H

#define ABS(num)  ((num) < 0.0 ? -(num) : (num))
#define MIN(a, b) (a < b ? a : b)

void imprimeMatriz(double *A, int tam);
void imprimeVetor(double *v, int tam);
void forwardSubstitution(double *A, double *x, double *b, int tam);
void retroSubstitution(double *A, double *x, double *b, int tam);
void trocaLinhas(double *A, double *b, int tam, int k, int l);
void metodoDeGauss(double *A, double *b, double *L, int tam);
void fatoracaoLU(double *A, double *L, int tam);

#endif