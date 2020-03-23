#include"stdio.h"
#include"omp.h"

int main()
{
	#pragma omp parallel
	{
		printf{"\n Hellow world !"};
	}
	printf{"\n program Exit !"};
}