#ifndef SLE_SOLVER
#define SLE_SOLVER

#include "matrix.hpp"
#include "vector.hpp"

template<class T>
matrix_t<T> join_mv(matrix_t<T>& A, vector_t<T>& b) {
	matrix_t<T>* M = new matrix_t<T>(A.get_rows, A.get_cols);
	
}

template<class T>
vector_t<T> sle_colver(matrix_t<T> A, vector_t<T> b) {
	
}

#endif // SLE_SOLVER

