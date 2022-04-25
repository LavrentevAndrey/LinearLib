#ifndef MATRIX_VECTOR_MULTIPLICATION
#define MATRIX_VECTOR_MULTIPLICATION

#include "matrix.hpp"
#include "vector.hpp"

template <class T>
vector_t<T> operator * (const matrix_t<T>& lM, const vector_t<T> &r) {
	if (lM.get_cols() != r.get_dim()) throw std::invalid_argument("Matrix_t and vector_t size multimplication error");
	int M = lM.get_rows(), N = lM.get_cols();
	T* tmp = new T[M]{};
	for (int i = 0; i < M; i++) {
		const T* a = lM.get_pointer() + i * N;
		for (int j = 0; j < N; j++) {
			tmp[i] += a[j] * r.get_cord(j);
		}
	}
	vector_t<T> res(tmp, M);
	delete[] tmp;
	return res;
}

template<class T>
vector_t<T> operator * (const vector_t<T> &l, const matrix_t<T> &rM) {
	if (rM.get_rows() != l.get_dim()) throw std::invalid_argument("Matrix_t and vector_t size multimplication error");
	int M = rM.get_rows(), N = rM.get_cols();
	T* tmp = new T[N]{};
	for (int i = 0; i < M; i++) {
		T cur = l.get_cord(i);
		const T* a = rM.get_pointer() + i * N;
		for (int j = 0; j < N; j++) {
			tmp[j] += a[j] * cur;
		}
	}
	vector_t<T> res(tmp, N);
	delete[] tmp;
	return res;
}

#endif