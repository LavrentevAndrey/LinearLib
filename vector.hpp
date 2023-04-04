#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <math.h>
#include <vector>

template<class T>
class vector_t {
public:
	// Default constructor
	vector_t();
	
	// Zero constructor
	explicit vector_t(int n);

	// Initializer list constructor
	vector_t(std::initializer_list<T> cords, int n);

	// Constructor by array
	vector_t(T* cords, int n);

	// Constructor by vector
	vector_t(std::vector<T> cords);

	// Copy cnstructor
	vector_t(vector_t& object);

	// Move constructor
	vector_t(vector_t&& object) noexcept; //noexcept - for calm type casting

	//Assignment copy constructor
	vector_t& operator = (vector_t& object);

	// Assignment transfer constructor
	vector_t& operator = (vector_t&& object) noexcept; //noexcept - for calm type casting

	// Destructor
	~vector_t();

	inline T get_cord(int n) const;

	inline int get_dim() const;

	inline void set_cord(T value, int n);

	inline double norm() const;

	// v /= |v|
	inline void normalise();

	inline vector_t<T> get_norm() const;


	// Operator == overlaod for l-values
	bool operator == (vector_t<T>& object);

	// Operator == overlaod for r-values
	bool operator == (vector_t<T>&& object);

	// Operator += overload for l-values
	vector_t& operator += (vector_t<T>& object);
	
	// Operator += overload for r-values
	vector_t& operator += (vector_t<T>&& object);

	// Operator -= overload for l-values
	vector_t& operator -= (vector_t<T>& object);

	// Operator -= overload for r-values
	vector_t& operator -= (vector_t<T>&& object);

	// Operator *= overload for l-values
	T operator *= (vector_t<T>& object);
	
	// Operator *= overload for r-values
	T operator *= (vector_t<T>&& object);

	// Operator *= overload for r-value coeficient
	vector_t& operator *= (const T ratio);

	// Operator /= overload for r-value coeficient
	vector_t& operator /= (const T ratio);

	// Vector sum overload
	template<class U> friend vector_t<U> operator + (vector_t<U> lhs, vector_t<U> rhs);

	// Vector minus overload
	template<class U> friend vector_t<U> operator - (vector_t<U> lhs, vector_t<U> rhs);

	// Dot product overload
	template<class U> friend U operator * (vector_t<U> lv, vector_t<U> rv);

	// Vector product overload
	template<class U> friend vector_t<U> operator % (vector_t<U> lhs, vector_t<U> rhs);

	// Out overload
	template <class U> friend std::ostream& operator<< (std::ostream& os, const vector_t<U>& M);

private:
	std::vector<T> data;
	int dim;
};

//------------------------------------------------------------------
//------------------------------------------------------------------

template<class T>
vector_t<T>::vector_t() {
	data = std::vector<T>();
	dim = 0;
}

template<class T>
vector_t<T>::vector_t(int n) {
	data = std::vector<T>(n);
	dim = n;
}

template<class T>
vector_t<T>::vector_t(std::initializer_list<T> cords, int n) {
	data = std::vector<T>(n);
	int i = 0;
	for (T elem : cords) {
		if (i >= n)
			break;
		data[i] = elem;
		i++;
	}
	dim = n;
}

template<class T>
vector_t<T>::vector_t(T* cords, int n) {
	data = std::vector<T>(n);
	for (int i = 0; i < n; i++)
		data[i] = cords[i];
	dim = n;
}

template<class T>
vector_t<T>::vector_t(std::vector<T> cords) {
	data = cords;
	dim = cords.size();
}

template<class T>
vector_t<T>::vector_t(vector_t& object) {
	data = object.data;
	dim = object.dim;
}

template<class T>
vector_t<T>::vector_t(vector_t&& object) noexcept { //noexcept - for calm type casting
	data = object.data;
	object.data = nullptr;
	dim = object.dim;
}

template<class T>
vector_t<T>& vector_t<T>::operator = (vector_t& object) {
	data = object.data;
	dim = object.dim;
	return *this;
}

template<class T>
vector_t<T>& vector_t<T>::operator = (vector_t&& object) noexcept { //noexcept - for calm type casting
	data = object.data;
	object.data = nullptr;
	dim = object.dim;
	return *this;
}

template<class T>
vector_t<T>::~vector_t() {}

//------------------------------------------------------------------
//------------------------------------------------------------------

template<class T>
bool vector_t<T>::operator == (vector_t<T>& object) {
	if (dim != object.dim)
		return false;

	for (int i = 0; i < object.dim; i++)
		if (fabs((double)(data[i] - object.data[i])) > 1e-9)
			return false;

	return true;
}

template<class T>
bool vector_t<T>::operator == (vector_t<T>&& object) {
	if (dim != object.dim)
		return false;

	for (int i = 0; i < object.dim; i++)
		if (fabs((double)(data[i] - object.data[i])) > 1e-9)
			return false;

	return true;
}

template<class T>
vector_t<T>& vector_t<T>::operator += (vector_t<T>& object) {
	if (dim != object.dim)
		throw std::invalid_argument("Dimensions need to be the same!");

	for (int i = 0; i < dim; i++) {
		data[i] += object.data[i];
	}
	return *this;
}

template<class T>
vector_t<T>& vector_t<T>::operator += (vector_t<T>&& object) {
	if (dim != object.dim)
		throw std::invalid_argument("Dimensions need to be the same!");

	for (int i = 0; i < dim; i++) {
		data[i] += object.data[i];
	}
	return *this;
}

template<class T>
vector_t<T>& vector_t<T>::operator -= (vector_t<T>& object) {
	if (dim != object.dim)
		throw std::invalid_argument("Dimensions need to be the same!");

	for (int i = 0; i < dim; i++) {
		data[i] -= object.data[i];
	}
	return *this;
}

template<class T>
vector_t<T>& vector_t<T>::operator -= (vector_t<T>&& object) {
	if (dim != object.dim)
		throw std::invalid_argument("Dimensions need to be the same!");

	for (int i = 0; i < dim; i++) {
		data[i] -= object.data[i];
	}
	return *this;
}

template<class T>
T vector_t<T>::operator *= (vector_t<T>& object) {
	if (dim != object.dim)
		throw std::invalid_argument("Dimensions need to be the same!");

	T sum = static_cast<T>(0.0);
	for (int i = 0; i < dim; i++) {
		sum += data[i] * object.data[i];
	}
	return sum;
}

template<class T>
T vector_t<T>::operator *= (vector_t<T>&& object) {
	if (dim != object.dim)
		throw std::invalid_argument("Dimensions need to be the same!");

	T sum = static_cast<T>(0.0);
	for (int i = 0; i < dim; i++) {
		sum += data[i] * object[i];
	}
	return sum;
}

template<class T>
vector_t<T>& vector_t<T>::operator *= (const T ratio) {
	for (int i = 0; i < dim; i++) {
		data[i] *= ratio;
	}
	return *this;
} 

template<class T>
vector_t<T>& vector_t<T>::operator /= (const T ratio) {
	if (ratio == static_cast<T>(0.0))
		throw std::invalid_argument("Zero division, floating point exception!");

	for (int i = 0; i < dim; i++) {
		data[i] /= ratio;
	}
	return *this;
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template<class T>
T vector_t<T>::get_cord(int n) const {
	return data.at(n);
} // Numeration starts with sero

template<class T>
int vector_t<T>::get_dim() const {
	return dim;
}

template<class T>
void vector_t<T>::set_cord(T value, int n) {
	data.at(n) = value;
} // Numeration starts with sero

template<class T>
double vector_t<T>::norm() const {
	double sum = 0.0;
	for (T elem : data)
		sum += static_cast<double>(elem * elem);
	return sqrt(sum);
}

template<class T>
void vector_t<T>::normalise() {
	if (this->norm() != 0.0)
		*this /= this->norm();
}

template<class T>
vector_t<T> vector_t<T>::get_norm() const {
	double div = static_cast<T>(this->norm());
	if (div == 0.0)
		return vector_t<T>(std::vector<T>(dim));
	std::vector<T> out = data;
	for (int i = 0; i < dim; i++)
		out[i] /= div;
	return vector_t<T>(out);
}

//------------------------------------------------------------------
//------------------------------------------------------------------

template<class T>
vector_t<T> operator + (vector_t<T> lv, vector_t<T> rv) {
	if (lv.dim != rv.dim)
		throw std::invalid_argument("Dimensions need to be the same!");

	std::vector<T> out;
	out.reserve(rv.dim + 1);

	for (int i = 0; i < rv.dim; i++) {
		out.push_back(lv.data[i] + rv.data[i]);
	}

	return vector_t<T>(out);
}

template<class T>
vector_t<T> operator - (vector_t<T> lv, vector_t<T> rv) {
	if (lv.dim != rv.dim)
		throw std::invalid_argument("Dimensions need to be the same!");

	std::vector<T> out;
	out.reserve(rv.dim + 1);

	for (int i = 0; i < rv.dim; i++) {
		out.push_back(lv.data[i] - rv.data[i]);
	}

	return vector_t<T>(out);
}

template<class T>
vector_t<T> operator * (vector_t<T> lv, const T value) {
	std::vector<T> out;
	for (T elem : lv.data) {
		out.push_back(elem * value);
	}
	return out;
}

template<class T>
T operator * (vector_t<T> lv, vector_t<T> rv) {
	if (lv.dim != rv.dim)
		throw std::invalid_argument("Dimensions need to be the same!");

	T sum = static_cast<T>(0.0);
	for (int i = 0; i < lv.dim; i++) {
		sum += lv.data[i] * rv.data[i];
	}
	return sum;
}

template<class T>
vector_t<T> operator % (vector_t<T> lhs, vector_t<T> rhs) {
	if (lhs.dim != rhs.dim)
		throw std::invalid_argument("Dimensions need to be the same!");
	
	if (lhs.dim == 3) {
		return vector_t<T>({ lhs.data[1] * rhs.data[2] - lhs.data[2] * rhs.data[1], \
							lhs.data[0] * rhs.data[2] - lhs.data[2] * rhs.data[0], \
							lhs.data[0] * rhs.data[1] - lhs.data[1] * rhs.data[0] }, 3);
	}
	else
		throw std::invalid_argument("Unavailable option, if you want cross product of dimensions different from 3, you need another function!");
}

template<class T>
std::ostream& operator << (std::ostream& os, const vector_t<T>& M) {
	for (T elem : M.data)
		os << elem << " ";
	os << std::endl;
	return os;
}

#endif