#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <chrono>
#include <cmath>

#include <lapacke.h>


// чтение матрицы из файла
int read_matrix(float *v, std::string str, int size)
{

    std::ifstream file(str);

    if (!file.is_open())
        return 0;


    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            file >> v[i * size + j];
        }
    }

    return 1;
}

// округление всех элементов матрицы
void matrix_round(float *v, int size)
{

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            float sign = std::copysignf(1, v[i * size + j]);
            v[i * size + j] = sign * round(std::copysignf(1, v[i * size + j]) * v[i * size + j] * 100) / 100;
        }
    }
}


// печать матрицы в файл 
void print_matrix_to_file(float *matrix, std::string name, int size)
{
    std::ofstream fout(name);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            fout << matrix[j * size + i] << " ";
        }
        fout << std::endl;
    }
    fout.close();
}



int main(int argc, char *argv[]) {
    
    
    int matrix_size = 256;
    std::string input_file_name = "matrix256.txt";
    std::string output_file_name = "tridiag_matrix256.txt";

    float *matrix = new float[matrix_size * matrix_size]; 
    if (!read_matrix(matrix, input_file_name, matrix_size))
    {
        std::cout << "Err open file";
        return 1;
    }


    float *D = new float[matrix_size];
    float *E = new float[matrix_size - 1];
    float *tau = new float[matrix_size];



    
    auto start = std::chrono::steady_clock::now();
    int info = LAPACKE_ssytrd(LAPACK_ROW_MAJOR, 'U', matrix_size, matrix, matrix_size, D, E, tau);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    matrix_round(matrix, matrix_size);
    print_matrix_to_file(matrix, output_file_name, matrix_size);

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    delete[] D;
    delete[] E;
    delete[] tau;
    delete[] matrix;


    return 0;
}