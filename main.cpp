#include "matrix2.hpp"

void gemm_v1(int M, int N, int K, const float* A, const float* B, float* C)
{
    for (int i = 0; i < M; ++i)
    {
        float* c = C + i * N;
        for (int j = 0; j < N; ++j)
            c[j] = 0;
        for (int k = 0; k < K; ++k)
        {
            const float* b = B + k * N;
            float a = A[i * K + k];
            for (int j = 0; j < N; ++j)
                c[j] += a * b[j];
        }
    }
}

int main() {
    float* A = new float[9];
    for (int i = 0; i < 9; i++)
        A[i] = rand() % 30;
    float* B = new float[9];
    for (int i = 0; i < 9; i++)
        B[i] = rand() % 120;
    float* C = new float[9];
    gemm_v1(2, 2, 2, A, B, C);
    Matrix2<float> Am(3, 3, A), Bm(3, 3, B), Cm(3, 3);
    std::cout << Am << std::endl << Bm << std::endl;
    Cm = Am * Bm;
    std::cout << "Matrix multiplication\n";
    std::cout << Cm << std::endl;
    std::cout << "Matrix multiplication by const\n";
    Cm = (Am * (float)3);
    std::cout << Cm << std::endl;
    std::cout << "Matrix addition\n";
    Cm = Am + Bm;
    std::cout << Cm << std::endl;
    std::cout << "Matrix subtraction\n";
    Cm = Bm -Am;
    std::cout << Cm << std::endl;
    std::cout << "Matrix equality\n";
    std::cout << (Am == Bm) << "   " << (Am == Am) << std::endl;
	return 0;
}

// Here will be tests :)

// 1 0  0 1  0 1
// 1 1  2 1  2 2