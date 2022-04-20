#include <iostream>
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
    float* A = new float[4];
    A[0] = 1, A[1] = 0, A[2] = 1, A[3] = 1;
    float* B = new float[4];
    B[0] = 0, B[1] = 1, B[2] = 2, B[3] = 1;
    float* C = new float[4];
    gemm_v1(2, 2, 2, A, B, C);
    Matrix2<float> Am(2, 2, A), Bm(2, 2, B), Cm(2, 2);
    Cm = Am * Bm;
    std::cout << Cm << std::endl;
    int i = 1;
	return 0;
}

// Here will be tests :)

// 1 0  0 1  0 1
// 1 1  2 1  2 2