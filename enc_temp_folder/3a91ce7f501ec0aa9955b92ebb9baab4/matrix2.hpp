#ifndef MATRIX2_H
#define MATRIX2_H

#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <math.h>
#include <vector>

template <class T>
class Matrix2 {
public:
	Matrix2();
	Matrix2(int rows, int cols);
	Matrix2(int rows, int cols, const T* data);
	Matrix2(const Matrix2<T>& matrix);
	Matrix2(int rows, int cols, const std::vector<T>* inputData);

	~Matrix2();

	// configuration methods
	bool resize(int rows, int cols);
	void set_to_identity();

	// get/set methods
	inline T get(int row, int col);
	inline bool set(int row, int col, const T value);
	int get_rows();
	int get_cols();

	// overloads of operators
	bool operator== (const Matrix2<T>& M);
	bool operator== (Matrix2<T>&& M);
	bool compare(const Matrix2<T>&& M, double tolerance);

	Matrix2<T>& operator = (const Matrix2<T>& M);
	
	template <class U> friend Matrix2<U> operator + (const Matrix2<U>& lM, const Matrix2<U>& rM);
	template <class U> friend Matrix2<U> operator + (const Matrix2<U>& lM, const U& r);
	template <class U> friend Matrix2<U> operator + (const U& l, const Matrix2<U>&rM);

	template <class U> friend Matrix2<U> operator - (const Matrix2<U>& lM, const Matrix2<U>& rM);
	template <class U> friend Matrix2<U> operator - (const Matrix2<U>& lM, const U& r);
	template <class U> friend Matrix2<U> operator - (const U& l, const Matrix2<U>& rM);

	template <class U> friend Matrix2<U> operator * (const Matrix2<U>& lM, const Matrix2<U>& rM);
	template <class U> friend Matrix2<U> operator * (const Matrix2<U>& lM, const U& r);
	template <class U> friend Matrix2<U> operator * (const U& l, const Matrix2<U>& rM);

	template <class U> friend std::ostream& operator<< (std::ostream& os, const Matrix2<U>& M);

	// computing inverse matrix
	bool inverse();

private:
	inline int pos_to_ind(int row, int col);
	inline bool is_square();
	inline bool close_enough(T n1, T n2);
	inline void mult_row(int i, T coef); // i = i * coef
	inline void mult_add_row(int i, int j, T coef); // i = i + j * coef
	inline void swap_row(int i, int j);
	bool separate(Matrix2<T>* M1, Matrix2<T>* M2, int col);
	bool join(const Matrix2<T>& M);
	inline int find_row_with_max_element(int col, int start_row);
	inline void print();

	T* m_data;
	int m_rows, m_cols, m_elem; 
};

//------------------------------------------------------------------
//------------------------------------------------------------------

// Default constructor
template<class T>
Matrix2<T>::Matrix2() {
	m_rows = 1;
	m_cols = 1;
	m_elem = 1;
	m_data = new T[m_elem];
	m_data[0] = 0.0;
}

// Empty MxN matrix
template<class T>
Matrix2<T>::Matrix2(int rows, int cols) {
	m_rows = rows;
	m_cols = cols;
	m_elem = rows * cols;
	m_data = new T[m_elem];
	for (int i = 0; i < m_elem; i++) {
		m_data[i] = 0.0;
	}
}

// Matrix from array
template<class T>
Matrix2<T>::Matrix2(int rows, int cols, const T* data) {
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
Matrix2<T>::Matrix2(const Matrix2<T>& matrix) {
	m_rows = matrix.m_rows;
	m_cols = matrix.m_cols;
	m_elem = matrix.m_elem;
	m_data = new T[m_elem];
	for (int i = 0; i < this->m_elem; i++) {
		this->m_data[i] = matrix.m_data[i];
	}
}

// IHAAA VECTORS!!!
template<class T>
Matrix2<T>::Matrix2(int rows, int cols, const std::vector<T>*  data) {
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
Matrix2<T>::~Matrix2() {
	if (m_data != nullptr)
		delete[] m_data;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

// Configuration
template<class T>
bool Matrix2<T>::resize(int rows, int cols) {
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
void Matrix2<T>::set_to_identity() {
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
T Matrix2<T>::get(int row, int col) {
	int ind = pos_to_ind(row, col);
	return (ind >= 0) ? this->m_data[ind] : throw std::invalid_argument("No such index in matrix");
}

// Put value in matrix by index
template<class T>
bool Matrix2<T>::set(int row, int col, const T value) {
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
int Matrix2<T>::get_rows() {
	return m_rows;
}

// This is rude too
template<class T>
int Matrix2<T>::get_cols() {
	return m_cols;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template<class T>
std::ostream& operator<< (std::ostream& os, Matrix2<T>& M) {
	int m = M.get_rows(), n = M.get_cols();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			os << M.get(i, j) << " ";
		os << std::endl;
	}
	os << std::endl;
	return os;
}

template<class T>
bool Matrix2<T>::operator== (const Matrix2<T>& M) {
	if (M.m_cols != this->m_cols || M.m_rows != this->m_rows)
		return false;

	for (int i = 0; i < this->m_elem; i++) {
		if (close_enough(this->m_data[i], M.m_data[i])) {
			return false;
		}
	}
	return true;
}

template<class T>
bool Matrix2<T>::operator== (Matrix2<T>&& M) {
	if (M.m_cols != this->m_cols || M.m_rows != this->m_rows)
		return false;
	
	for (int i = 0; i < this->m_elem; i++) {
		if (close_enough(this->m_data[i], M.m_data[i])) {
			return false;
		}
	}
	return true;
}

template<class T>
bool compare(const Matrix2<T> M, double coef) {
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
Matrix2<T>& Matrix2<T>::operator = (const Matrix2<T>& M) {
	this->m_cols = M.m_cols;
	this->m_rows = M.m_rows;
	this->m_elem = M.m_elem;
	this->m_data = new T[m_elem];
	for (int i = 0; i < m_elem; i++) {
		m_data[i] = M.m_data[i];
	}
	return *this;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template <class T>
Matrix2<T> operator+ (const Matrix2<T>& lM, const Matrix2<T>& rM) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM.m_data[i] + rM.m_data[i];

	Matrix2<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator+ (const Matrix2<T>& lM, const T& r) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM.m_data[i] + r;

	Matrix2<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator+ (const T& l, const Matrix2<T>& rM) {
	T* tmp = new T[rM.m_elem];
	for (int i = 0; i < rM.m_elem; i++)
		tmp[i] = rM.m_data[i] + l;

	Matrix2<T> res(rM.m_rows, rM.m_cols, tmp);
	delete[] tmp;
	return res;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template <class T>
Matrix2<T> operator- (const Matrix2<T>& lM, const Matrix2<T>& rM) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM.m_data[i] - rM.m_data[i];

	Matrix2<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator- (const Matrix2<T>& lM, const T& r) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM.m_data[i] - r;

	Matrix2<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator- (const T& l, const Matrix2<T>& rM) {
	T* tmp = new T[rM.m_elem];
	for (int i = 0; i < rM.m_elem; i++)
		tmp[i] = l - rM.m_data[i];

	Matrix2<T> res(rM.m_rows, rM.m_cols, tmp);
	delete[] tmp;
	return res;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template <class T>
Matrix2<T> operator* (const Matrix2<T>& lM, const T& r) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM.m_data[i] * r;

	Matrix2<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator* (const T& l, const Matrix2<T>& rM) {
	T* tmp = new T[rM.m_elem];
	for (int i = 0; i < rM.m_elem; i++)
		tmp[i] = l * rM.m_data[i];

	Matrix2<T> res(rM.m_rows, rM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator* (const Matrix2<T>& lM, const Matrix2<T>& rM) {
	// if (lM.m_cols != rM.m_rows) std::cout << "Matrix2 size multimplication error" << std::endl;
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
	Matrix2<T> res(M, N, tmp);
	delete[] tmp;
	return res;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

// Helping functions

// Separating two matrices with slicing them by a column number
template<class T>
bool Matrix2<T>::separate(Matrix2<T>* M1, Matrix2<T>* M2, int col) {
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
bool Matrix2<T>::join(const Matrix2<T>& M) {
	int cols1 = m_cols;
	int cols = this->m_cols + M.m_cols;
	
	if (M.m_rows != this->m_rows)
		throw std::invalid_argument("Attemp to join two unequal in rows matrices!");
	int rows = M.m_rows;
	

	this->print();

	T* M_data = new T[rows * cols];
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			int ind = i * cols + j;
			if (j < cols1)
				M_data[ind] = this->m_data[i * cols1 + j];
			else
				M_data[ind] = M.m_data[i * cols1 + j];
		}
	}

	Matrix2<T> A(3, 6, M_data);
	std::cout << A;

	m_elem = cols * rows;
	m_cols = cols;
	delete[] m_data;
	
	m_data = new T[m_elem];
	for (int i = 0; i < m_elem; i++)
		m_data[i] = M_data[i];
	//m_data = M_data;

	this->print();

	return true;
}

// Translate ij position in matrix into index of array
template<class T>
inline int Matrix2<T>::pos_to_ind(int row, int col) {
	return (row < m_rows && row >= 0 && col < m_cols && col >= 0) ? row * m_cols + col : -1;
}

template<class T>
inline bool Matrix2<T>::close_enough(T n1, T n2) {
	return fabs(n1 - n2) < 1e-9;
}

template<class T>
inline bool Matrix2<T>::is_square() {
	return m_cols == m_rows ? true : false;
}

template<class T>
inline void Matrix2<T>::swap_row(int i, int j) {
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
inline void Matrix2<T>::mult_add_row(int i, int j, T coef) {
	T* I = m_data + i * m_cols;
	T* J = m_data + j * m_cols;

	for (int k = 0; k < m_cols; k++) {
		I[k] += J[k] * coef;
	}
}

template<class T>
inline int Matrix2<T>::find_row_with_max_element(int col, int start_row) {
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
inline void Matrix2<T>::mult_row(int i, T coef) {
	T* I = m_data + i * m_cols;
	for (int k = 0; k < m_cols; k++)
		I[k] *= coef;
}

template<class T>
inline void Matrix2<T>::print() {
	for (int i = 0; i < m_rows; i++) {
		for (int j = 0; j < m_cols; j++)
			std::cout << this->get(i, j) << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template<class T>
bool Matrix2<T>::inverse() {
	if (!this->is_square())
		throw std::invalid_argument("Can't invert not square matrix!");
																												
	Matrix2<T> Ident(m_rows, m_cols);
	Ident.set_to_identity();																												
	int cols = m_cols;
	join(Ident);

	this->print();

	int cRow, cCol;
	int rowWithMaxElem;
	int count = 0;
	int maxIter = 100;
	bool flag = false, out = false;
	while (count < maxIter || flag == false) {
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

			for (int i = cRow + 1; i < m_rows; i++) {
				if (!close_enough(m_data[i * m_cols + cCol], 0.0)) {
					int rowOneIndex = cCol;

					T rowOneValue = m_data[rowOneIndex * m_cols + cCol];

					if (!close_enough(rowOneValue, 0.0)) {
						T correction = -(m_data[i * m_cols + cCol] / rowOneValue);
						this->mult_add_row(i, rowOneIndex, correction);
						//std::cout << "Multiply row: " << rowOneIndex << " by " << m_data[i * m_cols + cCol] << " and add to row " << i << std::endl;
					}
				}
			}

			for (int i = cCol + 1; i < cols; i++) {
				if (!close_enough(m_data[cRow * m_cols + i], 0.0)) {
					int rowOneIndex = i;

					T rowOneValue = m_data[rowOneIndex * m_cols + i];

					if (!close_enough(rowOneValue, 0.0)) {
						T correction = -(m_data[cRow * m_cols + i] / rowOneValue);
						this->mult_add_row(cRow, rowOneIndex, correction);
						//std::cout << "Multiply row: " << rowOneIndex << " by " << m_data[i * m_cols + cCol] << " and add to row " << cRow << std::endl;
					}
				}
			}

			flag = true;
			for (int i = 0; i < m_rows; i++) {
				for (int j = 0; j < cols; j++) {
					if (!(m_data[i * m_cols + j] == 0 && i != j || i == j && m_data[i * m_cols + j] == 1)) {
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
			}
		}
		count++;
	}
	return out;
}

#endif