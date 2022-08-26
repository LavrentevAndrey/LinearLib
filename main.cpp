#include <gtest/gtest.h>

#include "matrix.hpp"
#include "vector.hpp"
#include "marix_vector_multiplication.hpp"
#include "SLE_solver.hpp"

void test1() {
    double A[9], B[9], C[9];
    for (int i = 0; i < 9; i++)
        A[i] = (double)(rand() % 30);
    for (int i = 0; i < 9; i++)
        B[i] = (double)(rand() % 30);
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
    double A1[6], B1[6];
    for (int i = 0; i < 6; i++)
        A1[i] = (double)(rand() % 30);
    for (int i = 0; i < 6; i++)
        B1[i] = (double)(rand() % 30);
    Am = matrix_t<double>(3, 2, A1), Bm = matrix_t<double>(2, 3, B1);
    std::cout << Am << std::endl;
    std::cout << Bm << std::endl;
    Cm = Am * Bm;
    std::cout << "Not square matrix multiplication\n";
    std::cout << Cm << std::endl;
    Cm = Bm * Am;
    std::cout << "Not square matrix multiplication\n";
    std::cout << Cm << std::endl;
}

void test2() {
    double A [9];
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

    std::cout << s.get_norm() << std::endl;
    std::cout << s.norm() << std::endl;
    s.normalise();
    std::cout << s << std::endl;
    std::cout << s.norm() << std::endl;

}

void test4() {
    double A[6], B[6];
    for (int i = 0; i < 6; i++)
        A[i] = (double)(rand() % 10);
    matrix_t<double> Am = matrix_t<double>(2, 3, A);
    std::vector<double> vv1;
    for (int i = 0; i < 3; i++)
        vv1.push_back(rand() % 10);
    vector_t<double> v1(vv1);

    std::cout << Am << std::endl;
    Am.transpose();
    std::cout << Am << std::endl;
    Am.transpose();
    std::cout << Am << std::endl;
    std::cout << v1 << std::endl;

    std::cout << Am * v1;

    for (int i = 0; i < 6; i++)
        B[i] = (double)(rand() % 10);
    matrix_t<double> Bm = matrix_t<double>(3, 2, B);

    std::cout << "\n" << Bm << std::endl;
    std::cout << v1 << std::endl;
    std::cout << v1 * Bm;
}

void test5() {
    double B[5 * 5];
    for (int i = 0; i < 5 * 5; i++)
        B[i] = (double)(rand() % 10);
    B[10] = 17;
    matrix_t<double> Bm = matrix_t<double>(5, 5, B);
    matrix_t<double> Am = Bm;
    std::cout << "Entry matrix: \n" << Bm << std::endl;
    Bm.row_echelon();
    std::cout << "Entry matrix: \n" << Bm << std::endl;
    

    std::vector<double> vv1;
    for (int i = 0; i < 5; i++)
        vv1.push_back(rand() % 10);
    vector_t<double> v1(vv1);
    std::cout << "Entry vector: \n" << v1 << std::endl;

    std::vector<double> vv2 = { -1.001, -2.443, 0.61, 2.179, 0.831 };
    vector_t<double> v2(vv2);

    std::cout <<"Matrix x vector: \n" <<  Am * v2 << std::endl;

    vector_t<double> solution = sle_solver_NxN(Am, v1);
    std::cout << "Solution of matrix: \n" << solution << std::endl;

    double v[] = { 2, 1, -1, -1, 3, -1, 4, 2, -3 };
    
    matrix_t<double> M = matrix_t<double>(3, 3, v);
    //std::cout << M << std::endl;
    std::vector<double> vv3 = { 1, 2, 3};
    vector_t<double> v3(vv3);
    vector_t<double> solution1 = sle_solver_NxN(M, v3);
    std::cout << solution1 << std::endl;
    vector_t<double> sol({-1.0/7, 2.0/7, -1.0}, 3);
    assert(sol == solution1);
}

int main(int argc, char* argv[]) {
    // testing::InitGoogleTest(&argc, argv);
	// return RUN_ALL_TESTS();
    test5();
    return 0;
}

// Here will be tests :)

// 1 0  0 1  0 1
// 1 1  2 1  2 2
