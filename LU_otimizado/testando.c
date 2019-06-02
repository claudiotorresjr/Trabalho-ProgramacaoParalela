#include <stdio.h>

#include <omp.h>

int main(int argc, char const *argv[])
{	
	#pragma omp parallel num_threads(2)
	{	
		int ID = omp_get_thread_num();
		for (int i = ID; i < 10; i+=2)
		{
			printf("th %d pegando a linha %d\n", ID, i);
		}
	}
	return 0;
}