#include<iostream>
#include<vector>
#include <fstream>
#include <cmath>
#include <ctime> 
#include <chrono>
#include <stdexcept>
//#include <cblas.h>


// сюда кладу все векторы для дальнейшей проверки
std::vector<std::vector<float>> vectors_database;

// Функиция выполняет перемножение матриц PAP
int pap_transformation(std::vector<float>& a, std::vector<float>& b,
    std::vector<float>& matrix, float digit)
{
    if (a.size() != b.size() || a.size() != sqrt(matrix.size()))
        return 1;

    int size = a.size();

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrix[i * size + j] -= a[i] * b[j] * digit;
            matrix[i * size + j] -= a[j] * b[i] * digit;
        }
    }
    return 0;
}


// печать матрицы в файл 
void print_matrix_to_file(std::vector<float>& matrix, std::string name)
{
    int size = sqrt(matrix.size());
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

// печать вектора в файл
void print_vector_to_file(std::vector<float>& vector_name, std::string name)
{
    int size = vector_name.size();
    std::ofstream fout;
    fout.open("file.txt", std::ios::app);

    for (int i = 0; i < size; i++)
    {
        fout << vector_name[i] << " ";
    }
    fout << std::endl;
    fout.close();
}



// создание пустой матрицы
std::vector<float> create_matrix(int size)
{
    std::vector<float> a(size * size);
    return a;
}

// скалярное умножение векторов
int dot(std::vector<float>& a, std::vector<float>& b, float& ans)
{
    if (a.size() != b.size())
        return 1;

    ans = 0;
    int size = a.size();

    for (int i = 0; i < size; i++)
        ans += a[i] * b[i];
    return 0;
}

int dot_matrix_vector(std::vector<float>& matrix, std::vector<float>& b, std::vector<float>& ans)
{
    int size = ans.size();
    //Процедура xGEMM вычисляет следующее выражение: С = αAB + βC
    //CblasRowMajor - матрицы хранятся по строкам
    //CblasNoTrans - ничего не делать
    //size - размер матрицы

    //cblas_sgemv(CblasRowMajor, CblasNoTrans, size, size, 1, &matrix[0], size, &b[0], 1, 1.0f, &ans[0], 1);

    for (int k = 0; k < size; k++)
    {
        ans[k] = 0;
        for (int i = 0; i < size; i++)
            ans[k] += matrix[i + k * size] * b[i];
    }
    return 0;
}


// преобразование матрицы по алгоритму Хаусхолдера
int hausholders_transformation(std::vector<float>& matrix)
{
    int size = sqrt(matrix.size());
    std::ofstream outfile("log.txt");
    for (int i = 0; i < size - 2; i++)
    {
        float norm = 0;
        for (int j = i + 1; j < size; j++)
        {
            norm += matrix[size * j + i] * matrix[size * j + i];
        }
        norm = sqrt(norm);

        float s = std::copysignf(1, matrix[(i + 1) * size + i]) * norm;
        float r = sqrt(2 * matrix[(i + 1) * size + i] * s + 2 * s * s);

        std::vector<float> w(size);
        w[i + 1] = (matrix[(i + 1) * size + i] + s) * (1 / r);

        for (int j = i + 2; j < size; j++)
        {
            if (j >= i + 2)
            {
                w[j] = matrix[j * size + i] * (1 / r);
            }
            if (j <= i)
            {
                w[j] = 0;
            }
        }

        //the road to hell is paved with good intensions

        // сохранение вектора в файл log
        int size = w.size();
        for (int i = 0; i < size; i++)
        {
            outfile << w[i] << " ";
        }
        outfile << std::endl;

        // сохранение вектора в массив для проверки
        //vectors_database.push_back(w);

        std::vector<float> v(size);
        if (dot_matrix_vector(matrix, w, v))
        {
            std::cout << "matrix-vector multiplication error!";
            return 1;
        }

        float c = 0;
        if (dot(v, w, c))
        {
            std::cout << "vector-vector multiplication error!";
            return 1;
        }


        std::vector<float> q(size);
        for (int j = 0; j < size; j++)
        {
            q[j] = v[j] - c * w[j];
        }

        pap_transformation(w, q, matrix, 2);
    }
    outfile.close();
    return 0;
}

// чтение матрицы из файла
int read_matrix(std::vector<float>& v, std::string str)
{

    int size = sqrt(v.size());
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
void matrix_round(std::vector<float>& v)
{
    int size = sqrt(v.size());
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            float sign = std::copysignf(1, v[i * size + j]);
            v[i * size + j] = sign * round(std::copysignf(1, v[i * size + j]) * v[i * size + j] * 100) / 100;
        }
    }
}

int main()
{
    int matrix_size = 4;
    std::string input_file_name = "matrix4.txt";
    std::string output_file_name = "tridiag_matrix4.txt";

    auto matrix = create_matrix(matrix_size);
    if (!read_matrix(matrix, input_file_name))
    {
        std::cout << "Err open file";
        return 1;
    }
    auto start = std::chrono::steady_clock::now();

    hausholders_transformation(matrix);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    matrix_round(matrix);
    print_matrix_to_file(matrix, output_file_name);

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}

