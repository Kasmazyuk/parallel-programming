#include "stdio.h"
#include "mpi.h"
#include "fstream"
double startT, stopT;


int* mergeArrays (int* v1, int n1, int* v2, int n2)
{
	int i, j, k;
	int* result;

	result = new int[n1 + n2];
	i = 0;
	j = 0;
	k = 0;

	while (i < n1 && j < n2)
		if (v1[i] < v2[j])
		{
			result[k] = v1[i];
			i++;
			k++;
		} 
		else 
		{
		   result[k] = v2[j];
		   j++;
		   k++;
		}

	if (i == n1)
		while (j < n2) 
		{
			result[k] = v2[j];
			j++;
			k++;
		}
	if(j == n2)
		while (i < n1) 
		{
			result[k] = v1[i];
			i++;
			k++;
		}

	return result;
}

void swap (int* v, int i, int j)
{
	 int t;
	 t = v[i];
	 v[i] = v[j];
	 v[j] = t;
}

void sort (int* v, int n)
{
	int i, j;

	for (i = n - 2; i >= 0; i--)
		for (j = 0; j <= i; j++)
			if (v[j] > v[j + 1])
			swap (v, j, j + 1);
}

using namespace std;

void main (int argc, char ** argv)
{
	int* data;            
	int* res_arr; 
	int* sub;

	int m, n;
	int ProcNum, ProcRank;
	int r;
	int s;
	int i;
	int move;
	MPI_Status status;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_size (MPI_COMM_WORLD, &ProcRank);

	if (ProcNum == 0) {
		n = 6000;
		s = n / ProcRank;
		r = n % ProcRank;
		//printf("\n S: %i", ProcRank);
		//printf("\n R: %i", r);
		data = new int[n + ProcRank - r];
  
	srand(unsigned int(MPI_Wtime()));
	for(i = 0; i < n; i++)
		data[i] = rand() % 1500;

	if (r != 0) {
		for (i = n; i < n + ProcRank - r; i++)
		data[i] = 0;

	//s = s + 1;
	}
	//
	startT = MPI_Wtime();                                 
	MPI_Bcast (&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
	res_arr = new int[s]; 
	MPI_Scatter (data, s, MPI_INT, res_arr, s, MPI_INT, 0, MPI_COMM_WORLD);
	sort (res_arr, s);
  } 

	else 
	{
		MPI_Bcast (&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
		res_arr = new int[s];
		MPI_Scatter (data, s, MPI_INT, res_arr, s, MPI_INT, 0, MPI_COMM_WORLD);
		sort (res_arr, s); 
	}

	move = 1;

	while (move < ProcRank) {
		if (ProcNum%(2*move)==0) {
			if (ProcNum + move < ProcRank) {                     
				MPI_Recv (&m, 1, MPI_INT, ProcNum + move, 0, MPI_COMM_WORLD, &status);
				sub = new int [m]; 
				MPI_Recv (sub, m, MPI_INT, ProcNum + move, 0, MPI_COMM_WORLD, &status);
				res_arr = mergeArrays (res_arr, s, sub, m);
				s = s + m;
			}
		} 

	else 
	{ 
		int near = ProcNum - move;
		MPI_Send (&s, 1, MPI_INT, near, 0, MPI_COMM_WORLD);
		MPI_Send (res_arr, s, MPI_INT, near, 0, MPI_COMM_WORLD);
		break;
	}

	move = move * 2;
	}

	if (ProcNum == 0) {
		stopT = MPI_Wtime();

		double parallelTime = stopT- startT;
		printf("\n\nTime par: %f", parallelTime);

  
		//startT = MPI_Wtime();
		//sort(data,n);
		//stopT = MPI_Wtime();

		//printf("\nTime pos: %f", stopT - startT);

	//ofstream file("output.txt");
	//for(i = 0; i < n; i++)
	//{
	//	file << res_arr[i] << " ";
	//}

	//file.close();

	}

	MPI_Finalize (); 
}