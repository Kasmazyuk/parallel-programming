#include "mpi.h"
#include <cstdlib>
#include <string>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

int main (int argc, char **argv)
{
	srand(time(0));

	int n = 10;
	int borderL, borderR;
	char s = 'E';
	char *str;
	int par_count = 0;
	int count = 0;
	double timest, timefin;

	str = new char [n];
	for (int i = 0; i < n; i++)
	{
		str[i] = rand()%10 + 65;
	  //printf("%c", str[i]);
	}

	int ProcNum, ProcRank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Status Status;

	int k = n / ProcNum;
	int proc_count = 0;

	MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank == 0) 
		timest = MPI_Wtime();

	MPI_Bcast(str, n, MPI_CHAR, 0, MPI_COMM_WORLD);

	borderL = (int)(k*ProcRank);
	borderR = (int)(k*(ProcRank + 1));
	if (ProcRank == ProcNum - 1) 
		borderR = n;
	for (int j = borderL; j < borderR; j++)
		if (str[j] == s)
			proc_count++;
	MPI_Reduce(&proc_count, &par_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
	{
		timefin = MPI_Wtime();
		float res = (float)par_count/(float)n;
		printf("Repeatabillity = %f\n", res);
		printf("Time work = %f", timefin - timest);

	}

	MPI_Finalize();
	return 0;
}