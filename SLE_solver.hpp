#ifndef SLE_SOLVER
#define SLE_SOLVER

#include "matrix.hpp"
#include "vector.hpp"

template<class T>
matrix_t<T> join_mv(matrix_t<T>& A, vector_t<T>& b) {
	int M = A.get_rows(), N = A.get_cols();
	if (M != b.get_dim())
		throw std::invalid_argument("Not the same vector and matrix dimensions in join_mv");

	T* data = new T[M * (N + 1)];
	T* row;
	for (int i = 0; i < M; i++) {
		row = data + i * (N + 1);
		for (int j = 0; j < N; j++)
			row[j] = A.get(i, j);
		row[N] = b.get_cord(i);
	}

	matrix_t<T> out(M, N + 1, data);
	return out;
}

template<class T>
vector_t<T> sle_solver_NxN(matrix_t<T> A, vector_t<T> b) {
	int M = A.get_rows(), N = A.get_cols();
	if (M != N || b.get_dim() != M)
		throw std::invalid_argument("Cant solve, there are more variables than equations");

	matrix_t<T> partly_solved_sistem = (join_mv(A, b)).row_echelon();
	std::cout << partly_solved_sistem << std::endl;
	
	vector_t<T> out(M);
	T sum;
	for (int i = M - 1; i >= 0; i--) {
		sum = partly_solved_sistem.get(i, N);

		for (int j = N - 1; j > i; j--) {
			sum -= out.get_cord(j) * partly_solved_sistem.get(i, j);
		}

		out.set_cord(sum / partly_solved_sistem.get(i, i), i);
	}

	return out;
}

#endif // SLE_SOLVER

