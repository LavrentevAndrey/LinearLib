#ifndef MATRIX2_H
#define MATRIX2_H

template <class T>
class Matrix2 {
public:
	Matrix2();
	Matrix2(int rows, int cols);
	Matrix2(int rows, int cols, const T* data);
	Matrix2(const Matrix2<T>& matrix);

	~Matrix2();

	// configuration methods
	bool resize(int rows, int cols);

	// get/set methods
	T get(int row, int col);
	bool set(int row, int col, const T value);
	int get_rows();
	int get_cols();

	// overloads of operators
	bool operator== (const Matrix2<T>& M);
	
	template <class U> friend Matrix2<U> operator + (const Matrix2<U>& lM, const Matrix2<U>& rM);
	template <class U> friend Matrix2<U> operator + (const Matrix2<U>& lM, const U& r);
	template <class U> friend Matrix2<U> operator + (const U& l, const Matrix2<U>&rM);

	template <class U> friend Matrix2<U> operator - (const Matrix2<U>& lM, const Matrix2<U>& rM);
	template <class U> friend Matrix2<U> operator - (const Matrix2<U>& lM, const U& r);
	template <class U> friend Matrix2<U> operator - (const U& l, const Matrix2<U>& rM);

	template <class U> friend Matrix2<U> operator * (const Matrix2<U>& lM, const Matrix2<U>& rM);
	template <class U> friend Matrix2<U> operator * (const Matrix2<U>& lM, const U& r);
	template <class U> friend Matrix2<U> operator * (const U& l, const Matrix2<U>& rM);

private:
	int Pos2Ind(int row, int col);

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
	} // hmm the exception is asking for... (data size is lower then m_elem)
}

// Matrix from matrix
template<class T>
Matrix2<T>::Matrix2(const Matrix2<T>& matrix) {
	m_rows = matrix.m_rows;
	m_cols = matrix.m_cols;
	m_elem = matrix.m_elem;
	m_data = new T[m_elem];
	for (int i = 0; i < m_elem; i++) {
		m_data[i] = matrix.m_data[i];
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
	delete[] m_data;
	m_data = new T[m_elem];
	if (m_data != nullptr) {
		for (int i = 0; i < m_elem; i++)
			m_data[i] = 0.0;

		return true;
	}
	else
		return false;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

// Get value by position
template<class T>
T Matrix2<T>::get(int row, int col) {
	int ind = Pos2Ind(row, col);
	if (ind >= 0)
		return m_data[ind];
	else
		return 0;
}

// Put value in matrix by index
template<class T>
bool Matrix2<T>::set(int row, int col, const T value) {
	int ind = Pos2Ind(row, col);
	if (ind >= 0) {
		m_data[ind] = value;
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
bool Matrix2<T>::operator== (const Matrix2<T>& M) {
	if (M.m_cols != this->m_cols || M.m_rows != this->m_rows)
		return false;

	for (int i = 0; i < this->m_elem; i++) {
		if (this->m_data[i] != M.m_data[i]) {
			return false;
		}
	}
	return true;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template <class T>
Matrix2<T> operator+ (const Matrix2<T>& lM, const Matrix2<T>& rM) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM[i] + rM[i];

	Matrix2<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator+ (const Matrix2<T>& lM, const T& r) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM[i] + r;

	Matrix2<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator+ (const T& l, const Matrix2<T>& rM) {
	T* tmp = new T[rM.m_elem];
	for (int i = 0; i < rM.m_elem; i++)
		tmp[i] = rM[i] + l;

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
		tmp[i] = lM[i] - rM[i];

	Matrix2<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator- (const Matrix2<T>& lM, const T& r) {
	T* tmp = new T[lM.m_elem];
	for (int i = 0; i < lM.m_elem; i++)
		tmp[i] = lM[i] - r;

	Matrix2<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator- (const T& l, const Matrix2<T>& rM) {
	T* tmp = new T[rM.m_elem];
	for (int i = 0; i < rM.m_elem; i++)
		tmp[i] = l - rM[i];

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
		tmp[i] = lM[i] * r;

	Matrix2<T> res(lM.m_rows, lM.m_cols, tmp);
	delete[] tmp;
	return res;
}

template <class T>
Matrix2<T> operator* (const T& l, const Matrix2<T>& rM) {
	T* tmp = new T[rM.m_elem];
	for (int i = 0; i < rM.m_elem; i++)
		tmp[i] = l * rM[i];

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

// Translate ij position in matrix into index of array
template<class T>
int Matrix2<T>::Pos2Ind(int row, int col) {
	if (row < m_rows && row >= 0 && col < m_cols && col >= 0)
		return row * m_cols + col;
	else
		return -1;
}

#endif