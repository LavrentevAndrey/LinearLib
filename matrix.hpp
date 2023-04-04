#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <math.h>
#include <vector>
#include "vector.hpp"

template <class T>
class matrix_t {
public:
	matrix_t();
	matrix_t(int rows, int cols);
	matrix_t(int rows, int cols, const T* data);
	matrix_t(const matrix_t<T>& matrix);
	matrix_t(matrix_t<T>&& matrix);
	matrix_t(int rows, int cols, const std::vector<T>* inputData) noexcept;

	virtual ~matrix_t();

	// configuration methods
	bool resize(int rows, int cols);
	void set_to_identity();

	// get/set methods
	inline T get(int row, int col) const;
	inline bool set(int row, int col, const T value);
	inline int get_rows() const;
	inline int get_cols() const;
	inline T* get_pointer() const;

	// overloads of operators
	bool operator== (const matrix_t<T>& M);
	bool operator== (const matrix_t<T>&& M);
	bool compare(const matrix_t& M, double tolerance);

	const matrix_t& operator = (const matrix_t& M);
	const matrix_t& operator = (matrix_t<T>&& M);
	
	template <class U> friend matrix_t<U> operator + (const matrix_t<U>& lM, const matrix_t<U>& rM);
	template <class U> friend matrix_t<U> operator + (const matrix_t<U>& lM, const U& r);
	template <class U> friend matrix_t<U> operator + (const U& l, const matrix_t<U>&rM);

	template <class U> friend matrix_t<U> operator - (const matrix_t<U>& lM, const matrix_t<U>& rM);
	template <class U> friend matrix_t<U> operator - (const matrix_t<U>& lM, const U& r);
	template <class U> friend matrix_t<U> operator - (const U& l, const matrix_t<U>& rM);

	template <class U> friend matrix_t<U> operator * (const matrix_t<U>& lM, const matrix_t<U>& rM);
	template <class U> friend matrix_t<U> operator * (const matrix_t<U>& lM, const U& r);
	template <class U> friend matrix_t<U> operator * (const U& l, const matrix_t<U>& rM);


	template <class U> friend std::ostream& operator<< (std::ostream& os, const matrix_t<U>& M);

	// useful functions
	inline bool is_square();
	inline bool is_row_echelon();
	bool separate(matrix_t* M1, matrix_t* M2, int col);
	bool join(const matrix_t& M);

	// computing inverse matrix
	bool inverse();

	// Transposing matrix
	void transpose();

	// computing determinant
	T determinant(matrix_t M);
	T LU_determinant();

	// returning row echelon form of matrix
	void row_echelon();

	int get_rank() const;

private:
	inline int pos_to_ind(int row, int col) const;
	inline bool close_enough(T n1, T n2) const;
	inline void mult_row(int i, T coef); // i = i * coef
	inline void mult_add_row(int i, int j, T coef); // i = i + j * coef
	inline void swap_row(int i, int j);
	inline int find_row_with_max_element(int col, int start_row);
	inline void print();

	T* m_data;
	int m_rows, m_cols, m_elem; 
};

//------------------------------------------------------------------
//------------------------------------------------------------------

// Default constructor
template<class T>
matrix_t<T>::matrix_t() {
	m_rows = 1;
	m_cols = 1;
	m_elem = 1;
	m_data = new T[m_elem];
	m_data[0] = 0.0;
}

// Empty MxN matrix
template<class T>
matrix_t<T>::matrix_t(int rows, int cols) {
	m_rows = rows;
	m_cols = cols;
	m_elem = rows * cols;
	m_data = new T[m_elem];
	for (int i = 0; i < m_elem; i++) {
		m_data[i] = 0;
	}
}

// Matrix from array
template<class T>
matrix_t<T>::matrix_t(int rows, int cols, const T* data) {
	m_rows = rows;
	m_cols = cols;
	m_elem = rows * cols;
	m_data = new T[m_elem];
	for (int i = 0; i < m_elem; i++) {
		m_data[i] = data[i];
	}
}

// Matrix from matrix
template<class T>
matrix_t<T>::matrix_t(const matrix_t<T>& matrix) {
	m_rows = matrix.m_rows;
	m_cols = matrix.m_cols;
	m_elem = matrix.m_elem;
	m_data = new T[m_elem];
	for (int i = 0; i < this->m_elem; i++) {
		this->m_data[i] = matrix.m_data[i];
	}
}

template<class T>
matrix_t<T>::matrix_t(matrix_t<T>&& matrix) {
	m_rows = matrix.m_rows;
	m_cols = matrix.m_cols;
	m_elem = matrix.m_elem;
	m_data = matrix.m_data;
	matrix.m_data = nullptr;
}

// IHAAA VECTORS!!!
template<class T>
matrix_t<T>::matrix_t(int rows, int cols, const std::vector<T>*  data) noexcept {
	//assert(data->size() >= cols * rows); // hmm the exception is asking for... but is it so nessesary?
	m_rows = rows;
	m_cols = cols;
	m_elem = rows * cols;
	m_data = new T[m_elem];
	for (int i = 0; i < m_elem; i++) {
		m_data[i] = data->at(i);
	} 
}

// KILL 'EM ALL!!!
template<class T>
matrix_t<T>::~matrix_t() {
	if (m_data != nullptr)
		delete[] m_data;
	m_data = nullptr;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

// Configuration
template<class T>
bool matrix_t<T>::resize(int rows, int cols) {
	m_rows = rows;
	m_cols = cols;
	m_elem = rows * cols;
	delete[] this->m_data;
	this->m_data = new T[m_elem];
	if (this->m_data != nullptr) {
		for (int i = 0; i < m_elem; i++)
			this->m_data[i] = 0.0;

		return true;
	}
	else
		return false;
}
template<class T>
void matrix_t<T>::set_to_identity() {
	if (!this->is_square())
		throw std::invalid_argument("Cannot form an identity matrix that is not square!");
	for (int i = 0; i < m_rows; i++) {
		for (int j = 0; j < m_cols; j++) {
			if (i == j)
				this->m_data[i * m_cols + j] = 1.0;
			else
				this->m_data[i * m_cols + j] = 0.0;
		}
	}
}

//------------------------------------------------------------------
//------------------------------------------------------------------

// Get value by position
template<class T>
T matrix_t<T>::get(int row, int col) const {
	int ind = pos_to_ind(row, col);
	return (ind >= 0) ? this->m_data[ind] : throw std::invalid_argument("No such index in matrix");
}

// Put value in matrix by index
template<class T>
bool matrix_t<T>::set(int row, int col, const T value) {
	int ind = pos_to_ind(row, col);
	if (ind >= 0) {
		this->m_data[ind] = value;
		return true;
	}
	else
		return false;
}

// How rude, it is not permissible to deal in this form (Give me this, give me that)
template<class T>
int matrix_t<T>::get_rows() const {
	return m_rows;
}

// This is rude too
template<class T>
int matrix_t<T>::get_cols() const {
	return m_cols;
}

template<class T>
T* matrix_t<T>::get_pointer() const {
	return m_data;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template<class T>
std::ostream& operator<< (std::ostream& os, const matrix_t<T>& M) {
	int m = M.get_rows(), n = M.get_cols();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			os << std::fixed << std::setprecision(3) << std::setw(6) << M.get(i, j) << " ";
		os << std::endl;
	}
	os << std::endl;
	return os;
}

template<class T>
bool matrix_t<T>::operator== (const matrix_t<T>& M) {
	if (M.m_cols != this->m_cols || M.m_rows != this->m_rows)
		return false;

	for (int i = 0; i < this->m_elem; i++) {
		if (!close_enough(this->m_data[i], M.m_data[i])) {
			return false;
		}
	}
	return true;
}

template<class T>
bool matrix_t<T>::operator== (const matrix_t<T>&& M) {
	if (M.m_cols != this->m_cols || M.m_rows != this->m_rows)
		return false;
	
	for (int i = 0; i < this->m_elem; i++) {
		if (!close_enough(this->m_data[i], M.m_data[i])) {
			return false;
		}
	}
	return true;
}

template<class T>
bool matrix_t<T>::compare(const matrix_t<T> &M, double coef) {
	if (M.m_cols != this->m_cols || M.m_rows != this->m_rows)
		return false;

	double errorSum = 0;
	for (int i = 0; i < this->m_elem; i++) {
		T elem1 = M.m_data[i];
		T elem2 = this->m_data[i];
		errorSum += (elem1 - elem2) * (elem1 - elem2);
	}
	
	double error = sqrt(errorSum / ((this->m_rows * this->m_cols) + 1));
	return error < coef;
}

template<class T>
const matrix_t<T>& matrix_t<T>::operator = (const matrix_t<T>& M) {
	if (this == &M)
            return *this;
	m_cols = M.m_cols;
	m_rows = M.m_rows;
	m_elem = M.m_elem;
	m_data = new T[m_elem];
	for (int i = 0; i < m_elem; i++) {
		m_data[i] = M.m_data[i];
	}
	return *this;
}

template<class T>
const matrix_t<T>& matrix_t<T>::operator = (matrix_t<T>&& M) {
	if (this == &M)
            return *this;
	delete[] m_data;
	m_cols = M.m_cols;
	m_rows = M.m_rows;
	m_elem = M.m_elem;
	m_data = M.m_data;
	M.m_data = nullptr;
	return *this;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template <class T>
matrix_t<T> operator+ (const matrix_t<T>& lM, const matrix_t<T>& rM) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM.m_data[i] + rM.m_data[i];

	matrix_t<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
matrix_t<T> operator+ (const matrix_t<T>& lM, const T& r) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM.m_data[i] + r;

	matrix_t<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
matrix_t<T> operator+ (const T& l, const matrix_t<T>& rM) {
	T* tmp = new T[rM.m_elem];
	for (int i = 0; i < rM.m_elem; i++)
		tmp[i] = rM.m_data[i] + l;

	matrix_t<T> res(rM.m_rows, rM.m_cols, tmp);
	delete[] tmp;
	return res;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template <class T>
matrix_t<T> operator- (const matrix_t<T>& lM, const matrix_t<T>& rM) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM.m_data[i] - rM.m_data[i];

	matrix_t<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
matrix_t<T> operator- (const matrix_t<T>& lM, const T& r) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM.m_data[i] - r;

	matrix_t<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
matrix_t<T> operator- (const T& l, const matrix_t<T>& rM) {
	T* tmp = new T[rM.m_elem];
	for (int i = 0; i < rM.m_elem; i++)
		tmp[i] = l - rM.m_data[i];

	matrix_t<T> res(rM.m_rows, rM.m_cols, tmp);
	delete[] tmp;
	return res;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template <class T>
matrix_t<T> operator* (const matrix_t<T>& lM, const T& r) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM.m_data[i] * r;

	matrix_t<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
matrix_t<T> operator* (const T& l, const matrix_t<T>& rM) {
	T* tmp = new T[rM.m_elem];
	for (int i = 0; i < rM.m_elem; i++)
		tmp[i] = l * rM.m_data[i];

	matrix_t<T> res(rM.m_rows, rM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
matrix_t<T> operator* (const matrix_t<T>& lM, const matrix_t<T>& rM) {
	if (lM.m_cols != rM.m_rows) throw std::invalid_argument("Matrix_t and matrix_t size multimplication error");
	int M = lM.m_rows, N = rM.m_cols, K = lM.m_cols;
	T* tmp = new T[M * N];
	for (int i = 0; i < M; i++) {
		T* c = tmp + i * N;
		for (int j = 0; j < N; j++)
			c[j] = 0;
		for (int k = 0; k < K; k++) {
			const T* b = rM.m_data + k * N;
			T a = lM.m_data[i * K + k];
			for (int j = 0; j < N; j++)
				c[j] += a * b[j];
		}
	}
	matrix_t<T> res(M, N, tmp);
	delete[] tmp;
	return res;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

// Helping functions

// Separating two matrices with slicing them by a column number
template<class T>
bool matrix_t<T>::separate(matrix_t<T>* M1, matrix_t<T>* M2, int col) {
	int rows = this->m_rows;
	int cols = this->m_rows;
	if (col <= 0 || col >= cols)
		return false;
	M1->resize(rows, col);
	M2->resize(rows, cols - col);

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (i < col)
				M1->set(i, j, get(i,j));
			else
				M2->set(i, j - col, get(i, j));
		}
	}
	return true;
}

// Joining two matrices in one by columns
template<class T>
bool matrix_t<T>::join(const matrix_t<T>& M) {
	int cols1 = m_cols;
	int cols = this->m_cols + M.m_cols;
	
	if (M.m_rows != this->m_rows)
		throw std::invalid_argument("Attemp to join two unequal in rows matrices!");
	int rows = M.m_rows;

	T* M_data = new T[rows * cols];
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int ind = i * cols + j;
			if (j < cols1)
				M_data[ind] = this->m_data[i * cols1 + j];
			else
				M_data[ind] = M.m_data[i * cols1 + j - 3];
		}
	}

	m_elem = cols * rows;
	m_cols = cols;
	delete[] m_data;
	m_data = M_data;

	return true;
}

// Translate ij position in matrix into index of array
template<class T>
inline int matrix_t<T>::pos_to_ind(int row, int col) const {
	return (row < m_rows && row >= 0 && col < m_cols && col >= 0) ? row * m_cols + col : -1;
}

template<class T>
inline bool matrix_t<T>::close_enough(T n1, T n2) const {
	return fabs(n1 - n2) < 1e-9;
}

template<class T>
inline bool matrix_t<T>::is_square() {
	return m_cols == m_rows ? true : false;
}

template <class T>
bool matrix_t<T>::is_row_echelon() {
	int lol;
	for (int i = 1; i < m_rows; ++i) {
		lol = i * m_cols;
		for (int j = 0; j < i; ++j) {
			if (!close_enough(m_data[lol + j], 0.0))
				return false;
		}
	}
	return true;
}

template<class T>
inline void matrix_t<T>::swap_row(int i, int j) {
	T tmp;
	T* I = m_data + i * m_cols;
	T* J = m_data + j * m_cols;

	for (int k = 0; k < m_cols; k++) {
		tmp = I[k];
		I[k] = J[k];
		J[k] = tmp;
	}
}

template<class T>
inline void matrix_t<T>::mult_add_row(int i, int j, T coef) {
	T* I = m_data + i * m_cols;
	T* J = m_data + j * m_cols;

	for (int k = 0; k < m_cols; k++) {
		I[k] += J[k] * coef;
	}
}

template<class T>
inline int matrix_t<T>::find_row_with_max_element(int col, int start_row) {
	T max = m_data[0];
	int row = 0, ind;
	
	for (int k = start_row; k < m_rows; k++) {
		ind = k * m_cols + col;
		if (fabs(m_data[ind]) > fabs(max)) {
			max = m_data[ind];
			row = k;
		}
	}

	return row;
}

template<class T>
inline void matrix_t<T>::mult_row(int i, T coef) {
	T* I = m_data + i * m_cols;
	for (int k = 0; k < m_cols; k++)
		I[k] *= coef;
}

template<class T>
inline void matrix_t<T>::print() {
	for (int i = 0; i < m_rows; i++) {
		for (int j = 0; j < m_cols; j++)
			std::cout << std::fixed << std::setprecision(3) << std::setw(6) << this->get(i, j) << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template<class T>
void matrix_t<T>::transpose() {
	if (m_rows ==  m_cols) {
		T tmp;
		for (int i = 1; i < m_rows; i++) {
			for (int j = 0; j < i; j++) {
				tmp = m_data[i*m_cols + j];
				m_data[i*m_cols + j] = m_data[j*m_rows + i];
				m_data[j*m_rows + i] = tmp;
			}
		}
	}
	else {
		int tmpi = m_rows;
		T* buf = new T[m_rows * m_cols];
		for (int i = 0; i < m_rows; i++) {
			for (int j = 0; j < m_cols; j++) {
				buf[j*m_rows + i] = m_data[i*m_cols + j];
			}
		}
		delete[] m_data;
		m_data = buf;
		m_rows = m_cols;
		m_cols = tmpi;
	}
}

template<class T>
bool matrix_t<T>::inverse() {
	if (!this->is_square())
		throw std::invalid_argument("Can't invert not square matrix!");
																												
	matrix_t Ident(m_rows, m_cols);
	Ident.set_to_identity();																												
	int cols = m_cols;
	join(Ident);

	//this->print();

	int cRow, cCol;
	int rowWithMaxElem;
	int count = 0;
	int maxIter = m_rows * 10;
	bool flag = false, out = false;
	while (count < maxIter && flag == false) {
		for (int diagIndex = 0; diagIndex < m_rows; diagIndex++) {
			cRow = diagIndex;
			cCol = diagIndex;

			rowWithMaxElem = this->find_row_with_max_element(cCol, cRow);
			if (rowWithMaxElem != cRow) {
				this->swap_row(rowWithMaxElem, cRow);
				//std::cout << "Swap rows: " << cRow << " and " << rowWithMaxElem << std::endl;
			}

			T cValue = m_data[cRow * m_cols + cCol];
			if (cValue != 1.0)
				this->mult_row(cRow, (T)1 / cValue);

			//this->print();

			for (int i = cRow + 1; i < m_rows; i++) {
				if (!close_enough(m_data[i * m_cols + cCol], 0.0)) {
					int rowOneIndex = cCol;

					T rowOneValue = m_data[rowOneIndex * m_cols + cCol];

					if (!close_enough(rowOneValue, 0.0)) {
						T correction = -(m_data[i * m_cols + cCol] / rowOneValue);

						//std::cout << "Multiply row: " << rowOneIndex << " by " << m_data[i * m_cols + cCol] << " and add to row " << i << std::endl;

						this->mult_add_row(i, rowOneIndex, correction);

						//this->print();
					}
				}
			}

			//this->print();

			for (int i = cCol + 1; i < cols; i++) {
				if (!close_enough(m_data[cRow * m_cols + i], 0.0)) {
					int rowOneIndex = i;

					T rowOneValue = m_data[rowOneIndex * m_cols + i];

					if (!close_enough(rowOneValue, 0.0)) {
						T correction = -(m_data[cRow * m_cols + i] / rowOneValue);

						//std::cout << "Multiply row: " << rowOneIndex << " by " << m_data[i * m_cols + cCol] << " and add to row " << cRow << std::endl;

						this->mult_add_row(cRow, rowOneIndex, correction);

						//this->print();
					}
				}
			}

			//this->print();

			flag = true;
			for (int i = 0; i < m_rows; i++) {
				for (int j = 0; j < cols; j++) {
					if ((!close_enough(m_data[i * m_cols + j], 0.0) || i == j) && (i != j || m_data[i * m_cols + j] != 1)) {
						flag = false;
					}
				}
			}

			if (flag == true) {
				out = true;
				T* L = new T[m_rows * cols];
				for (int i = 0; i < m_rows; i++)
					for (int j = 0; j < cols; j++)
						L[i * cols + j] = m_data[i * m_cols + cols + j];

				m_cols = cols;
				m_elem = cols * m_rows;
				delete[] m_data;
				m_data = new T[m_elem];
				for (int i = 0; i < m_elem; i++)
					m_data[i] = L[i];
				delete[] L;
				break;
			}
		}
		count++;
	}
	return out;
}

template<class T>
T matrix_t<T>::determinant(const matrix_t<T> M) {
	if (M.m_cols == 2)
		return M.m_data[0] * M.m_data[3] - M.m_data[1] * M.m_data[2];
	int cols = M.m_cols, rowSmal, row;
	T sum = 0;
	T* tmp;
	tmp = new T[(cols - 1) * (cols - 1)];
	
	for (int i = 0; i < cols; i++) {

		// submatrix array
		for (int j = 1; j < cols; j++) {
			rowSmal = (j - 1) * (cols - 1);
			row = j * cols;
			for (int k = 0; k < i; k++)
				tmp[rowSmal + k] = M.m_data[row + k];
			for (int k = i + 1; k < cols; k++) {
				tmp[rowSmal + k - 1] = M.m_data[row + k];
			}
		}

		matrix_t subMatrix(cols - 1, cols - 1, tmp);
		if (i % 2 == 0) 
			sum += M.m_data[i] * determinant(subMatrix);
		else
			sum -= M.m_data[i] * determinant(subMatrix);
	}

	delete[] tmp;

	if (cols == 1)
		return M.m_data[0];
	return sum;
}

template<class T>
T matrix_t<T>::LU_determinant() {
	if (!this->is_square())
		throw std::invalid_argument("Can't invert not square matrix!");

	int n = m_cols;
	std::vector<std::vector<T>> L(n, std::vector<double>(n)), U(n, std::vector<double>(n));

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			U[i][j] = m_data[i * n + j];

	for (int i = 0; i < n; i++)
		for (int j = i; j < n; j++)
			L[j][i] = U[j][i] / U[i][i];

	for (int k = 1; k < n; k++)
	{
		for (int i = k - 1; i < n; i++)
			for (int j = i; j < n; j++)
				L[j][i] = U[j][i] / U[i][i];

		for (int i = k; i < n; i++)
			for (int j = k - 1; j < n; j++)
				U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
	}

	T det = (T)1.0;
	for (int i = 0; i < n; i++)
		det *= U[i][i];
	return det;
}

template <class T>
void matrix_t<T>::row_echelon() {
	if (m_cols < m_rows)
		throw std::invalid_argument("The matrix must have at least as many columns as rows.");

	int cRow, cCol;
	int maxCount = m_cols	;
	int count = 0;
	bool completeFlag = false;
	while ((!completeFlag) && (count < maxCount)) {
		for (int diagIndex = 0; diagIndex < m_rows; diagIndex++) {
			cRow = diagIndex;
			cCol = diagIndex;

			// std::cout << "Row: " << m_data[diagIndex * m_cols + diagIndex] << " " << diagIndex << std::endl;
			// this->print();
			if (m_data[diagIndex * m_cols + diagIndex] == 0) {
				for (int i = cRow + 1; i < m_rows; i++) {
					if (!close_enough(m_data[i * m_cols + cCol], 0)) {
						// std::cout << "Swap rows: " << cRow << " and " << i << std::endl;
						this->swap_row(i, cRow);
						// this->print();
						break;
					}
				}
			}

			for (int i = cRow + 1; i < m_rows; ++i) {
				if (!close_enough(m_data[pos_to_ind(i, cCol)], 0.0)) {
					int rowOneIndex = cCol;
					T rowOneValue = m_data[pos_to_ind(rowOneIndex, cCol)];
					if (!close_enough(rowOneValue, 0.0)) {
						T correctionFactor = -(m_data[pos_to_ind(i, cCol)] / rowOneValue);
						mult_add_row(i, rowOneIndex, correctionFactor);
					}
					//std::cout << "Step " << i << std::endl;
					//this->print();
				}
			}
		}
		completeFlag = this->is_row_echelon();
		count++;
		// std::cout << count << std::endl;
	}
	// this->print();
}

template<class T>
int matrix_t<T>::get_rank() const {
	int rank = std::max(m_cols,m_rows), size = m_cols * m_rows;
	std::vector<bool> line_used (m_rows);
	T *data = new T[size];
	for (int i = 0; i < size; i++)
		data[i] = m_data[i];

	for (int i = m_cols - 1; i >= 0; i--) {
		int j;
		for (j = 0; j < m_rows; j++)
			if (!line_used[j] && !close_enough(data[j*m_cols + i], 0.0))
				break;
		if (j == m_rows)
			rank--;
		else {
			line_used[j] = true;
			for (int p = i - 1; p >= 0; p--)
			 	data[j*m_cols + p] /= data[j*m_cols + i];
			for (int k = m_rows - 1; k >= 0; k--)
				if (k != j && !close_enough(data[k*m_cols + i], 0.0))
					for (int p = i - 1; p >= 0; p--)
						data[k*m_cols + p] -= data[j*m_cols + p] * data[k*m_cols + i];
		}
	}
	delete[] data;
	return rank;
}

#endif