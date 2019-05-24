#ifndef _MATRIZ_H
#define _MATRIZ_H

#define BLOCK_SIZE 8
#define Unroll 8

void escreveMatriz(double *C, unsigned int tam);

void multMatrizNormal(double *A, double *B, double *C, unsigned int tam);
void multMatrizTransposta(double *A, double *Bt, double *C, unsigned int tam);

void multMatrizNormalBloco(double *A, double *B, double *C, unsigned int tam);
void multMatrizTranspostaBloco(double *A, double *Bt, double *C, unsigned int tam);

#endif