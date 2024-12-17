
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
#define COUNT 1000000

// Функция для численного интегрирования по формуле прямоугольников
double rectangularIntegration(double a, double b, int n) {
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

    int a = A;
    int b = B;

    double res=0;
    double local_res=0;

    //Инициализаиця MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcessCount);

    int NumProcessCount = ProcessCount - 1;

    steady_clock::time_point start = high_resolution_clock::now();


    //Root процесс
    if (Rank == 0) {

        //Если процесс единственный
        if (ProcessCount == 1) {

            steady_clock::time_point start = high_resolution_clock::now();

            for (int i = 0; i < COUNT; i++) {
                rectangularIntegration(A, B, NumProcessCount - 1);
            }

            steady_clock::time_point end = high_resolution_clock::now();

            cout << "Single-process run\n" << duration_cast<milliseconds>(end - start).count() << " miliseconds" << endl;

            //Завершение MPI (в случае одного процесса)
            MPI_Finalize();

            return 0;

        }


    }

    for (int i = 0; i < COUNT/NumProcessCount; i++) {

        //Широковещательная отправка левой границы интервала
        MPI_Bcast(&a, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //Широковещательная отправка правой границы интервала
        MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //Широковещательная отправка количества numproc-процессов
        MPI_Bcast(&NumProcessCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
        //Работа внутри numproc-процесса
        if (Rank > 0) {
            // Длина подинтервала
            double len = ((b - a) / (double)NumProcessCount);
            // Левая граница подинтервала
            double local_a = a + (Rank - 1) * len;
            // Правая граница подинтервала
            double local_b = local_a + len;
            // Значение функции в подинтервале
            local_res = (exp(local_a) + exp(local_b)) / 2;
            local_res *= len;
        }

        //Слив значений из процессов в сумму
        MPI_Reduce(&local_res, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        //Барьерная синхронизация процессов
        MPI_Barrier(MPI_COMM_WORLD);

    }

    //Подсчет времени работы
    if (Rank == 0) {
        
        steady_clock::time_point end = high_resolution_clock::now();
        cout << "Multi-process run\n" << duration_cast<milliseconds>(end - start).count() <<" miliseconds" << endl;
        cout << "Amount of processes " << NumProcessCount << endl;
    }

    //Завершение MPI
    MPI_Finalize();
	
	return 0;
}