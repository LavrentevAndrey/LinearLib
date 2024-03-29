#ifndef SLE_SOLVER
#define SLE_SOLVER

#include "matrix.hpp"
#include "vector.hpp"

constexpr int INFINITE_SOLVES = -1;
constexpr int NO_SOLVES = -2;

template<class T>
matrix_t<T> join_mv(const matrix_t<T>& A, const vector_t<T>& b) {
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
	delete[] data;
	return out;
}

template<class T>
matrix_t<T> join_mv(const vector_t<T>& b, const matrix_t<T>& A) {
	int M = A.get_rows(), N = A.get_cols();
	std::cout << M << " iam in new\n";
	if (M != b.get_dim())
		throw std::invalid_argument("Not the same vector and matrix dimensions in join_mv");

	T* data = new T[M * (N + 1)];
	T* row;
	for (int i = 0; i < M; i++) {
		row = data + i * (N + 1);
		row[0] = b.get_cord(i);
		for (int j = 0; j < N; j++)
			row[j + 1] = A.get(i, j);
	}

	matrix_t<T> out(M, N + 1, data);
	delete[] data;
	return out;
}

template<class T>
int sle_solver_NxN(const matrix_t<T> &A, const vector_t<T> &b, vector_t<T> &result) {
	int M = A.get_rows(), N = A.get_cols();
	if (M != N || b.get_dim() != M)
		throw std::invalid_argument("Cant solve, there are more variables than equations");
	matrix_t<T> partly_solved_sistem(join_mv(A, b));
	partly_solved_sistem.row_echelon();
	std::cout << std::endl << partly_solved_sistem << std::endl;
	int big_rank = partly_solved_sistem.get_rank(), small_rank = A.get_rank();
	if ((big_rank == small_rank) && (small_rank < M)) {
		return INFINITE_SOLVES; 
	} else if (big_rank != small_rank) {
		return NO_SOLVES;
	}
	// std::cout << std::endl << partly_solved_sistem << std::endl;
	
	vector_t<T> out(M);
	T sum;
	if (partly_solved_sistem.get(0,0) != 0) {
		for (int i = M - 1; i >= 0; i--) {
			sum = partly_solved_sistem.get(i, N);

			for (int j = N - 1; j > i; j--) {
				sum -= out.get_cord(j) * partly_solved_sistem.get(i, j);
			}

			out.set_cord(sum / partly_solved_sistem.get(i, i), i);
		}
	}
	else if (partly_solved_sistem.get(M - 1, 0) != 0) {
		for (int i = 0; i < M; i++) {
			for (int j = N - 1; j >= N - i - 1; j--) {
				// std::cout << i << " " << j << std::endl;
				sum -= out.get_cord(N - j - 1) * partly_solved_sistem.get(i, j);
			}

			out.set_cord(sum / partly_solved_sistem.get(i, i), i);
		}
	}
	
	result = out;
	return 1;
}

#endif // SLE_SOLVER

