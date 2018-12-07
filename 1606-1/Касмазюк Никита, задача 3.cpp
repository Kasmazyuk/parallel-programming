#include <iostream>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <fstream>

using namespace std;


void PrintMatrix (double* matrix,int N){
	for (int i = 0; i < N; i++){
		for(int j = 0; j < N; j++)
			cout << matrix[i*N+j] << " ";
		cout << endl;
	}
	cout << endl;
}

void fPrintMatrix(double* matrix, int N, char name[20]){
	ofstream fout(name, ios_base::trunc);
	for (int i = 0; i < N; i++){
		for(int j = 0; j < N; j++)
			fout << matrix[i*N+j] << " ";
		fout << endl;
	}
	fout.close();
}

void RandMatrix(double* matrix1, double* matrix2, int N){
	srand (time (0));
	for (int i = 0; i < N; i++)
		for(int j = 0; j < N; j++){
			matrix1[i*N+j] = (double)rand() / 1000.0;
			matrix2[i*N+j] = (double)rand() / 1000.0;
		}
}

double* Trivial_alghorithm(double* matrix1, double* matrix2, int N){
	double* Rez = new double[N*N];
	double sum;
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++){
			sum = 0;
			for(int k = 0; k < N; k++)
				sum += matrix1[i*N+k] * matrix2[k*N+j];
			Rez[i*N+j]= sum;
		}
	return Rez;
}

double* Add(double* matrix1, double* matrix2, int N){
	double* Rez = new double[N*N];
	
	for (int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			Rez[i*N+j] = matrix1[i*N+j] + matrix2[i*N+j];

	return Rez;
}

double* Add(double* matrix1, double* matrix2, double* matrix3, double* matrix4, int N){
	double* Rez = new double[N*N];
	
	for (int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			Rez[i*N+j] = matrix1[i*N+j] + matrix2[i*N+j] + matrix3[i*N+j] + matrix4[i*N+j];

	return Rez;
}

double* Sub(double* matrix1, double* matrix2, int N){
	double* Rez = new double[N*N];
	
	for (int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			Rez[i*N+j] = matrix1[i*N+j] - matrix2[i*N+j];

	return Rez;
}

double* Sub(double* matrix1, double* matrix2, double* matrix3, double* matrix4, int N){//ñóììà 3-õ ìàòðèö âû÷åñòü ÷åòâåðòóþ
	double* Rez = new double[N*N];

	for (int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			Rez[i*N+j] = matrix1[i*N+j] + matrix2[i*N+j] + matrix3[i*N+j] - matrix4[i*N+j];

	return Rez;
}

double* alg(double* matrix1, double* matrix2, int N, int threshold){
	double* Rez;

	if(N <= threshold){
		Rez = Trivial_alghorithm(matrix1, matrix2, N);
	}
	else{
		Rez = new double[N*N];
		N = N/2;
		double* A[4];
		double* B[4];
		double* C[4];
		double* P[7];


		for(int i = 0; i < 4; i++){
			A[i] = new double[N*N];
			B[i] = new double[N*N];
		}


		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++){
				int index_new = i*N+j,index_old = 2*i*N+j, N_N = 2*N*N;
				A[0][index_new] = matrix1[index_old];
				A[1][index_new] = matrix1[index_old+N];
				A[2][index_new] = matrix1[index_old+N_N];
				A[3][index_new] = matrix1[index_old+N_N+N];

				B[0][index_new] = matrix2[index_old];
				B[1][index_new] = matrix2[index_old+N];
				B[2][index_new] = matrix2[index_old+N_N];
				B[3][index_new] = matrix2[index_old+N_N+N];
			}

		double *TMP = Add(A[0], A[3], N);
		double *_TMP = Add(B[0], B[3], N);
		P[0] = alg(TMP, _TMP, N,threshold);
		delete[] TMP;
		delete[] _TMP;

		TMP = Add(A[2], A[3], N);
		P[1] = alg(TMP, B[0], N, threshold);
		delete[] TMP;

		TMP = Sub(B[1], B[3], N);
		P[2] = alg(A[0], TMP, N, threshold);
		delete[] TMP;

		TMP = Sub(B[2], B[0], N);
		P[3] = alg(A[3], TMP, N, threshold);
		delete[] TMP;

		TMP = Add(A[0], A[1], N);
		P[4] = alg(TMP, B[3], N, threshold);
		delete[] TMP;

		TMP = Sub(A[2], A[0], N);
		_TMP = Add(B[0], B[1], N);
		P[5] = alg(TMP, _TMP, N, threshold);
		delete[] TMP;
		delete[] _TMP;

		TMP = Sub(A[1], A[3], N);
		_TMP = Add(B[2], B[3], N);
		P[6] = alg(TMP, _TMP, N, threshold);
		delete[] TMP;
		delete[] _TMP;


		C[0] = Sub(P[0], P[3], P[6], P[4], N);
		C[1] = Add(P[2], P[4], N);
		C[2] = Add(P[1], P[3], N);
		C[3] = Sub(P[0], P[2], P[5], P[1], N);


		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++){
					Rez[i*2*N+j] = C[0][i*N+j];
					Rez[i*2*N+j+N] = C[1][i*N+j];
					Rez[i*2*N+j+2*N*N] = C[2][i*N+j];
					Rez[i*2*N+j+2*N*N+N] = C[3][i*N+j];
			}

		for(int i = 0; i < 4; i++){
			delete[] A[i];
			delete[] B[i];
			delete[] C[i];
		}
		for(int i = 0; i < 7; i++)
			delete[] P[i];

	}

	return Rez;
}

int main(int argc, char** argv){
	double *matr_A, *matr_B, *matr_Rez_Seq, *matr_Rez_Par, 
		**A, **B, **TMP_Rez;
	int  N, thr = 64, sqr, new_N;

	int ProcSize = 1, ProcRank = 0;
	MPI_Status Status;

	double EndSeqAlg = 0; 
	double StartSeqAlg = 0;	
	double EndParAlg = 0; 
	double StartParAlg = 0;	 
	double TimeSeqAlg = 0;
	double TimeParAlg = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0){
		
		int k = 0;
		if (argc > 1)
			k = atoi (argv[1]);
		N = (int)pow(2.0,k);
		if(argc > 2)
			thr = atoi(argv[2]);

		
		matr_A= new double[N*N];
		matr_B= new double[N*N];
		matr_Rez_Par= new double[N*N];
		RandMatrix(matr_A, matr_B, N);


		if (N < 20){
			cout << "Matr A:" << endl;
			PrintMatrix(matr_A, N);
			cout << "Matr B:" << endl;
			PrintMatrix(matr_B, N);
		}
	

			// àëãîðèòì
		StartParAlg = MPI_Wtime();

		sqr = (int)sqrt((double)ProcSize), new_N = N/sqr;
		A = new double*[ProcSize], B = new double*[ProcSize];
	

		for(int i = 0; i < ProcSize; i++){
			A[i] = new double[N*N];
			B[i] = new double[N*N];
		}

		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++){
				A[sqr*(i/new_N)+j/new_N][(i%new_N)*new_N+(j%new_N)] = matr_A[i*N+j];
				B[sqr*(i/new_N)+j/new_N][(i%new_N)*new_N+(j%new_N)] = matr_B[i*N+j];
			}

		MPI_Bcast(&N, 1,  MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(&thr, 1,  MPI_INT, 0, MPI_COMM_WORLD);

		for(int i = 1; i < ProcSize; i++){
			int coef_A = sqr*(i / sqr), coef_B = i % sqr;
			for(int j = 0; j < sqr; j++){
				MPI_Send(A[coef_A], new_N*new_N,  MPI_DOUBLE, i , 0, MPI_COMM_WORLD); 
				MPI_Send(B[coef_B], new_N*new_N,  MPI_DOUBLE, i , 0, MPI_COMM_WORLD); 
				coef_A++;
				coef_B += sqr;
			}
		}

		for(int i = 0; i < sqr; i++){
			double* TMP = B[i];
			B[i] = B[i*sqr];
			B[i*sqr] = TMP;
		}
	}

	if(ProcRank !=0){
		MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		MPI_Bcast(&thr, 1, MPI_INT, 0, MPI_COMM_WORLD);

		sqr = (int)sqrt((double)ProcSize), new_N = N/sqr;

		A = new double*[sqr], B = new double*[sqr];
		for(int i = 0; i < sqr; i++){
			A[i] = new double[new_N*new_N];
			B[i] = new double[new_N*new_N];
		}

		for(int i = 0; i < sqr; i++){
			MPI_Recv(A[i], new_N*new_N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);
			MPI_Recv(B[i], new_N*new_N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);
		}
	}

	TMP_Rez = new double*[sqr+1];
	for(int i = 0; i < sqr; i++){
		TMP_Rez[i+1] = alg(A[i], B[i], new_N, thr);
	}
	if(ProcSize == 4)
		TMP_Rez[0] = Add(TMP_Rez[1], TMP_Rez[2], new_N);
	if(ProcSize == 16)
		TMP_Rez[0] = Add(TMP_Rez[1], TMP_Rez[2], TMP_Rez[3], TMP_Rez[4], new_N);

	if(ProcRank == 0){
		for(int i = 0; i < ProcSize; i++){
			delete[] A[i];
			delete[] B[i];
		}
	}
	else{
		for(int i = 0; i < sqr; i++){
			delete[] A[i];
			delete[] B[i];
		}
	}
	for(int i = 1; i < sqr+1; i++){
		delete[] TMP_Rez[i];
	}
	delete[] A;
	delete[] B;

	if (ProcRank != 0)	{
		MPI_Send(TMP_Rez[0], new_N*new_N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		delete[] TMP_Rez[0];
	}
	
	if (ProcRank == 0){
		int coef = sqrt(ProcSize);

		for(int i = 0; i < new_N; i++)
			for(int j = 0; j < new_N; j++)
				matr_Rez_Par[coef*i*new_N+j] = TMP_Rez[0][i*new_N+j];
		for(int k = 1; k < ProcSize; k++){

			MPI_Recv(TMP_Rez[0], new_N*new_N, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &Status);
			for(int i = 0; i < new_N; i++)
				for(int j = 0; j < new_N; j++)
					matr_Rez_Par[(k/coef)*new_N*N+(k%coef)*new_N+coef*i*new_N+j] = TMP_Rez[0][i*new_N+j];
		}

		//âûâîäû

		EndParAlg = MPI_Wtime();
		TimeParAlg = EndParAlg - StartParAlg;
		
		cout << "ParAlg time = " << TimeParAlg << endl;
		char fParAlg[20] = "ParAlg.txt";
			if (N < 20){
				cout << "Matr C:" << endl;
				PrintMatrix(matr_Rez_Par, N);
			}
			else{
				fPrintMatrix(matr_Rez_Par, N, fParAlg);
			}

		delete[] TMP_Rez[0];

		StartSeqAlg = MPI_Wtime();
		matr_Rez_Seq = alg(matr_A, matr_B, N, thr);
		EndSeqAlg = MPI_Wtime();
		TimeSeqAlg = EndSeqAlg - StartSeqAlg;
		
		cout << "SeqAlg time = " << TimeSeqAlg << endl;
		char fSeqAlg[20] = "SeqAlg.txt";
		if (N < 20){
			cout << "Matr C:" << endl;
			PrintMatrix(matr_Rez_Seq, N);
		}
		else{
			fPrintMatrix(matr_Rez_Seq, N, fSeqAlg);
		}

		//bool flag = false;
		//for(int k = 0; k < N; k++)
		//	for(int l = 0; l < N; l++)
		//		if(matr_Rez_Seq[k*N+l] != matr_Rez_Par[k*N+l]) 
		//			flag = true;


		//if(flag)
		//	cout << "matr_Rez_Seq != matr_Rez_Par" << endl;
		//if(!flag)
		//	cout << "matr_Rez_Seq == matr_Rez_Par" << endl;
		delete[] matr_Rez_Par;
		delete[] matr_Rez_Seq;
		delete[] matr_A;
		delete[] matr_B;
	}

	delete[] TMP_Rez;

	MPI_Finalize();

	return 0;
}
