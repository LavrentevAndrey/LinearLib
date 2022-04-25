#include "matrix.hpp"
#include "vector.hpp"
#include "marix_vector_multiplication.hpp"

void test1() {
    double* A = new double[9];
    for (int i = 0; i < 9; i++)
        A[i] = (double)(rand() % 30);
    double* B = new double[9];
    for (int i = 0; i < 9; i++)
        B[i] = (double)(rand() % 30);
    double* C = new double[9];
    matrix_t<double> Am(3, 3, A), Bm(3, 3, B), Cm(3, 3);
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
    Cm = Bm - Am;
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
    Am = matrix_t<double>(3, 2, A), Bm = matrix_t<double>(2, 3, B);
    std::cout << Am << std::endl;
    std::cout << Bm << std::endl;
    Cm = Am * Bm;
    std::cout << "Not square matrix multiplication\n";
    std::cout << Cm << std::endl;
}

void test2() {
    double* A = new double[9];
    for (int i = 0; i < 9; i++)
        A[i] = (double)(rand() % 30);
    matrix_t<double> Am(3, 3, A);

    std::cout << Am << std::endl << "Correct answer: " << Am.LU_determinant() << std::endl;
    std::cout << "Rec answ: " << Am.determinant(Am) << std::endl;
}

void test3() {
    std::vector<double> vv1, vv2;
    for (int i = 0; i < 3; i++)
        vv1.push_back(rand() % 10), vv2.push_back(rand() % 30);
    vector_t<double> v1(vv1), v2(vv2);
    std::cout << v1 << v2;
    vector_t<double> s(3);
    s = v1 + v2;
    std::cout << s << std::endl;
    s = v1 - v2;
    std::cout << s << std::endl;
    s = v1 % v2;
    std::cout << s << std::endl;
    std::cout << v1 * v2 << std::endl;
    s += v1;
    std::cout << s << std::endl;
    s -= v1;
    std::cout << s << std::endl;
    s *= 2;
    std::cout << s << std::endl;
    std::cout << (s *= v1) << std::endl;
    s /= 3;
    std::cout << s << std::endl;
    std::cout << (s == s) << std::endl;
    std::cout << (s == v1) << std::endl;
}

void test4() {
    double* A = new double[6];
    for (int i = 0; i < 6; i++)
        A[i] = (double)(rand() % 10);
    matrix_t<double> Am = matrix_t<double>(2, 3, A);
    std::vector<double> vv1;
    for (int i = 0; i < 3; i++)
        vv1.push_back(rand() % 10);
    vector_t<double> v1(vv1);

    std::cout << Am << std::endl;
    std::cout << v1 << std::endl;

    std::cout << Am * v1;

    double* B = new double[6];
    for (int i = 0; i < 6; i++)
        B[i] = (double)(rand() % 10);
    matrix_t<double> Bm = matrix_t<double>(3, 2, B);

    std::cout << "\n" << Bm << std::endl;
    std::cout << v1 << std::endl;
    std::cout << v1 * Bm;
}

int main() {
    test4();
	return 0;
}

// Here will be tests :)

// 1 0  0 1  0 1
// 1 1  2 1  2 2
