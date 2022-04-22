#include "matrix2.hpp"

int main() {
    double* A = new double[9];
    for (int i = 0; i < 9; i++)
        A[i] = (double)(rand() % 30);
    double* B = new double[9];
    for (int i = 0; i < 9; i++)
        B[i] = (double)(rand() % 30);
    double* C = new double[9];
    Matrix2<double> Am(3, 3, A), Bm(3, 3, B), Cm(3, 3);
    std::cout << Am << std::endl;
    Am.inverse();
    std::cout << Am << std::endl;
    std::cout << Bm << std::endl;
    Cm = Am * Bm;
    std::cout << "Matrix multiplication\n";
    std::cout << Cm << std::endl;
    std::cout << "Matrix multiplication by const\n";
    Cm = (Am * (double)3);
    std::cout << Cm << std::endl;
    std::cout << "Matrix addition\n";
    Cm = Am + Bm;
    std::cout << Cm << std::endl;
    std::cout << "Matrix subtraction\n";
    Cm = Bm -Am;
    std::cout << Cm << std::endl;
    std::cout << "Matrix equality\n";
    std::cout << (Am == Bm) << "   " << (Am == Am) << std::endl;
    Am.resize(3, 2);
    Bm.resize(2, 3);
    delete[] A, B, C;
    A = new double[6];
    for (int i = 0; i < 6; i++)
        A[i] = (double)(rand() % 30);
    B = new double[6];
    for (int i = 0; i < 6; i++)
        B[i] = (double)(rand() % 30);
    Am = Matrix2<double>(3, 2, A), Bm = Matrix2<double>(2, 3, B);
    std::cout << Am << std::endl;
    std::cout << Bm << std::endl;
    Cm = Am * Bm;
    std::cout << "Not square matrix multiplication\n";
    std::cout << Cm << std::endl;
    
	return 0;
}

// Here will be tests :)

// 1 0  0 1  0 1
// 1 1  2 1  2 2