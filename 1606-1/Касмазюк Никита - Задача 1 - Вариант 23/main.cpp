#include "mpi.h"
#include <cstdlib>
#include <string>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sstream>

using namespace std;

int main (int argc, char **argv)
{
	srand(time(0));

	int n;
	int borderL, borderR;
	char s;
	char *str;
	char *rbuf;
	int par_count = 0;
	int count = 0;
	double timest, timefin;

	stringstream convert(argv[1]); 
	stringstream convert1(argv[2]);

	if (!(convert >> n)) 
		n = 0; 

	if (!(convert1 >> s)) 
		s = '0'; 


	str = new char [n];
	rbuf = new char [n];
	//for (int i = 0; i < n; i++)
	//{
	//	str[i] = rand()%10 + 65;
	//  //printf("%c", str[i]);
	//}

	int ProcNum, ProcRank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Status Status;

	int g = n % ProcNum;
	int k = n / ProcNum;
	k = k + g;
	int proc_count = 0;

	//MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank == 0) 
	{
			/*str = new char [n];*/
			for (int i = 0; i < n; i++)
			{
				str[i] = rand()%10 + 65;
			  printf("%c ", str[i]);
	  
			}
			  printf("\n");
		/*for (int count=0; count < argc; ++count)
        std::cout << count << " " << argv[count] << '\n';*/

		timest = MPI_Wtime();
	}

	MPI_Scatter(str, k, MPI_CHAR, rbuf, k, MPI_CHAR, 0, MPI_COMM_WORLD);
	//MPI_Bcast(str, n, MPI_CHAR, 0, MPI_COMM_WORLD);

	//borderL = (int)(k*ProcRank);
	//borderR = (int)(k*(ProcRank + 1));
	//if (ProcRank == ProcNum - 1) 
	//	borderR = n;
	for (int j = 0; j < k; j++)
		if (rbuf[j] == s)
			proc_count++;
	MPI_Reduce(&proc_count, &par_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
	{
		timefin = MPI_Wtime();
		float res = (float)par_count/(float)n;
		printf("Repeatabillity = %f\n", res);
		printf("Time work      = %f", timefin - timest);

	}

	MPI_Finalize();
	return 0;
}