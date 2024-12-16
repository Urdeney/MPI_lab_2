
#include <iostream>
#include <iomanip>
#include <thread>
#include <chrono>
#include <cmath>
#include "mpi.h"

using namespace std;
using namespace std::chrono;

#define A 0
#define B 1


// Функция для численного интегрирования по формуле прямоугольников
double rectangularIntegrate(double a, double b, int n) {
    double len = (b - a) / n;
    double result = 0.0;

    for (int i = 0; i < n; ++i) {
        double local_a = a + i * len;
        double local_b = local_a + len;
        result += (exp(local_a) + exp(local_b))/2;
    }

    result *= (b-a) / n;
    return result;
}


int main()
{
    setlocale(LC_ALL, "Russian");
    srand(static_cast<unsigned int>(time(NULL)));
    
    int Rank = 0;
    int ProcessCount = 0;

    MPI_Status Status;

    int a = A;
    int b = B;

    double res;


    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcessCount);

    int WorkerProcessCount = ProcessCount - 1;

    steady_clock::time_point start = high_resolution_clock::now();


    //Root процесс
    if (Rank == 0) {

        //Если процесс единственный
        if (ProcessCount == 1) {

            rectangularIntegrate(A, B, 10);

            steady_clock::time_point end = high_resolution_clock::now();

            cout << "Single-process run\n" << duration_cast<milliseconds>(end - start).count() << endl;

            MPI_Finalize();

            return 0;

        }
    }


    //Отправка A
    MPI_Bcast(&a, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //Отправка B
    MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //Отправка N
    MPI_Bcast(&WorkerProcessCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //cout << "Hello i am working proccess" << endl;

    double len = ((b - a) / (double) WorkerProcessCount);
    double local_a = a + Rank * len;
    double local_b = local_a + len;
    double local_res= (exp(local_a) + exp(local_b)) / 2;
    local_res*= len;
    if (Rank == WorkerProcessCount) { local_res = 0;} //TODO убрать процесс- корень (лишний поток)

    
    cout << a << "\t" << b << "\t" << WorkerProcessCount << "\t" << Rank << endl;
    /*cout << "LOCAL len = " << len << endl;
    cout << "LOCAL A = " <<  local_a << endl;
    cout << "LOCAL B = " << local_b << endl;
    */
    cout << "LOCAL RES = " << local_res << endl;

    MPI_Reduce(&local_res, &res, 1, MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    
    if (Rank == 0) {
        cout <<"FINAL RES = " << res << endl;
        cout << rectangularIntegrate(A, B, 3) << endl;
    }

    MPI_Finalize();
	
	return 0;
}