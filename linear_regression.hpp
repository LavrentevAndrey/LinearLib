#ifndef LINEAR_REGRESSION_H
#define LINEAR_REGRESSION_H

#include "matrix.hpp"
#include "vector.hpp"
#include "SLE_solver.hpp"

template<class T>
int linear_regression(const matrix_t<T> &X_values, const vector_t<T> &y, vector_t<T> &result) {
    std::vector<T> data(X_values.get_rows(), (T)1.0);
    vector_t<T> ones(data);
    // std::cout << "Die here!\n";
    matrix_t<T> X(join_mv(ones, X_values));
    matrix_t<T> XT(X);
    XT.transpose();
    matrix_t<T> XTX(XT * X);
    // std::cout << "X: " << X << "XT: " << XT << "XTX: " << XTX << std::endl;
    // std::cout << "Y: " << y << std::endl;
    vector_t<T> XTY(XT * y);
    // std::cout << "XTY: " << XTY << std::endl;
    // std::cout << sle_solver_NxN(XTX, XTY, result) << std::endl;
    // std::cout << "That is it:\n" << result << std::endl;
    return sle_solver_NxN(XTX, XTY, result);
}

#endif // LINEAR_REGRESSION_H