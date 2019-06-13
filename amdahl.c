#include <stdio.h>

int main(int argc, char const *argv[])
{
	int proc;

	scanf("%d", &proc);
	
	printf("%lf\n", 1/(0.05 + (0.95/proc)));
	return 0;
}