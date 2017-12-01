#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <omp.h>

#define  Max(a,b) ((a)>(b)?(a):(b))

double maxeps = 0.1e-7;
int    itmax = 100;
double eps;


void relax(int size, double *matrix, int thread_num);
void init(int size, double *matrix, int thread_num);
void verify(int size, double *matrix,int thread_num); 

int 
main(int an, char **as)
{
	if (an < 3)
	{
		printf("Bad arguments, the number of arguments must be 2:\n\t1)size of matrix\n\t2)num of threads\n");
		return -1;
	}

	int size = atoi(as[1]);
    int thread_num = atoi(as[2]);
    double *A;
	int it;

    A = malloc(size * size * sizeof(double));

	double time = omp_get_wtime();
	init(size, A, thread_num);

	for(it = 1; it <= 100; it++)
	{
		eps = 0.0;
		relax(size, A, thread_num);
		//printf("Iteration = %d, Eps = %f\n", it, eps);
		if (eps < maxeps) break;
	}
	
	verify(size, A, thread_num);

	printf("Time %f\n\n", omp_get_wtime()-time);
	free(A);
	return 0;
}


void 
init(int size, double *matrix, int thread_num)
{ 
	int i,j;
#pragma omp parallel for shared(size, matrix, thread_num) private(i, j) num_threads(thread_num)
	for(i = 0; i <= size - 1; i++)
	{
		for(j = 0; j <= size - 1; j++)
		{
			if(i == 0 || i == size - 1 || j == 0 || j == size - 1)
			{
				matrix[i * size + j]= 0.;	
			}
			else
			{ 
				matrix[i * size + j]= ( 1. + i + j ) ;
			}
		}
	}
} 


void 
relax(int size, double *matrix, int thread_num)
{
	int i,j;
//#pragma omp declare reduction(maximum : double : omp_out = omp_in > omp_out ? omp_in : omp_out)
#pragma omp parallel for shared(size, matrix, thread_num) private(i,j) num_threads(thread_num) reduction(+:eps)
	for(i = 1; i <= size - 2; i++)
	{
		for(j = 1; j <= size-2; j++)
		{ 
			double e, tmp, dif;
			e = matrix[i * size + j];
			matrix[i * size + j] = (matrix[(i-1) * size + j] + matrix[(i+1) * size + j] +
									matrix[i * size + (j-1)] + matrix[i * size + (j+1)])/4.;
			eps = Max(eps, fabs(e-matrix[i * size + j]));
			tmp = fabs(e-matrix[i * size + j]);
			if ((dif = tmp - eps) > 0)   
			{
			 	eps += dif;
			}
		}
	}    
}


void 
verify(int size, double *matrix, int thread_num)
{ 
	int i,j;
	double s = 0.0;
#pragma omp parallel for reduction(+:s) shared(size, matrix, thread_num) private(i,j) num_threads(thread_num)
	for(i = 0; i <= size - 1; i++)
	{
		for(j = 0; j <= size-1; j++)
		{
			s += matrix[i * size + j] * (i+1) * (j+1) / (size * size);
		}
	}
	//printf("S = %f\n",s);
}


