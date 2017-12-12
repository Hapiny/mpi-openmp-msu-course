#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define  Max(a,b) ((a)>(b)?(a):(b))

double maxeps = 0.1e-7;
int    itmax = 100;
double eps;
double sum;

void relax(int num_rows, int begin_row);
void init(int num_rows, int begin_row);
void verify(int num_rows, int begin_row); 

#define SIZE 1024
int i, j;
int len_req, req_index; 

double (*A)[SIZE], (*B)[SIZE];
MPI_Request req[4];
MPI_Status st[4];
int myrank, num_rank;

void init(int num_rows, int begin_row);
void relax();
void verify();

int main(int an, char **as)
{
	
	int it;

	int begin_row, end_row, num_rows; 

	MPI_Init(&an, &as);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_rank);
	MPI_Barrier(MPI_COMM_WORLD);

	begin_row = (myrank * SIZE)/ num_rank;
	end_row = (((myrank + 1) * SIZE) / num_rank) - 1;
	num_rows = end_row - begin_row + 1;


	A = malloc((num_rows + 2) * SIZE * sizeof(double));
	B = malloc((num_rows) * SIZE * sizeof(double));

	double t1;
	if (myrank == 0) 
	{
		t1 = MPI_Wtime();
	}

	MPI_Barrier(MPI_COMM_WORLD);

	init(num_rows, begin_row);

	for(it = 1; it <= itmax; it++)
	{
		eps = 0.;
		relax(num_rows, begin_row);
		if (eps < maxeps) break;
	}

	verify(num_rows, begin_row);


	if (myrank == 0)
	{
		printf("Seconds = %f\n\n", MPI_Wtime()-t1);
	}
	MPI_Finalize();
	return 0;
}


void init(int num_rows, int begin_row)
{ 

	for(i = 1; i <= num_rows; i++)
	{
		for(j = 0; j <= SIZE-1; j++)
		{
			A[i][j] = 0.0;

			if((begin_row + i - 1 == 0) || (begin_row + i - 1 == SIZE-1) || (j == 0) || (j == SIZE-1))
				B[i-1][j]= 0.0;
			else 
				B[i-1][j]= (1. + begin_row + (i - 1) + j) ;
		}
	}
} 


void relax(int num_rows, int begin_row)
{

	for(i = 1; i <= num_rows; i++)
	{
		if(((i == 1) && (myrank == 0)) || ((i == num_rows) && (myrank == num_rank-1)))
			continue;
		for(j=1; j<= SIZE-2; j++)
		{
			A[i][j] = B[i-1][j];
		}
	}
	if (num_rank != 1)
	{
		if(myrank != 0)
			MPI_Irecv(&A[0][0], SIZE, MPI_DOUBLE, myrank - 1, 100, MPI_COMM_WORLD, &req[0]);
		if(myrank != num_rank - 1)
			MPI_Isend(&A[num_rows][0], SIZE, MPI_DOUBLE, myrank + 1, 100, MPI_COMM_WORLD, &req[2]);
		if(myrank != num_rank - 1)
			MPI_Irecv(&A[num_rows+1][0], SIZE, MPI_DOUBLE, myrank + 1, 101, MPI_COMM_WORLD, &req[3]);
		if(myrank != 0)
			MPI_Isend(&A[1][0], SIZE, MPI_DOUBLE, myrank - 1, 101, MPI_COMM_WORLD, &req[1]);
		
		len_req = 4; req_index = 0;
		if(myrank == 0)
		{
			len_req = 2;
			req_index = 2;
		}
		if(myrank == num_rank-1)
		{
			len_req = 2;
		}
		MPI_Waitall(len_req, &req[req_index], &st[0]);
	}

	double e;

	for(i = 1; i <= num_rows; i++)
	{
			if (((i==1) && (myrank == 0)) || ((i == num_rows) && (myrank == num_rank -1))) 
				continue;
			for(j = 1; j<= SIZE - 2; j++)
			{
				e = B[i-1][j];
				B[i-1][j] = (A[i-1][j] + A[i+1][j] + A[i][j-1] + A[i][j+1]) / 4.0;
				eps = Max(eps, fabs(e - B[i-1][j]));
			}
	}
}


void verify(int num_rows, int begin_row)
{ 
	double s=0.0;
	for(i = 0; i <= num_rows - 1; i++)
	{
		for(j = 0; j <= SIZE - 1; j++)
		{
			s = s + B[i][j] * (begin_row + (i - 1) + 1) * (j+1) / (SIZE * SIZE);
		}
	}
	MPI_Reduce(&s, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}





