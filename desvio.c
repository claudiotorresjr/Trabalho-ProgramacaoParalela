#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char const *argv[])
{
	FILE *arq;
	int i = 0;

	double vetor[20];

	double media = 0.0, sigma, p;
	int ch;

    arq = fopen("soma.tmp", "r");

    while(!feof(arq))
    {
    	fscanf(arq, "%lf\n", &vetor[i]);
    	media += vetor[i];
    	i++;
    }

    media /= 20.0;
    for(i = 0; i < 20; ++i)
    {
        p = p + pow(vetor[i] - media, 2);
    }
    sigma = sqrt(p/(19));
    printf("%.5f %.5f\n", media, sigma);

    return 0;
}