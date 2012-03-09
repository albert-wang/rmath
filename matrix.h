#ifndef MATRIX_HPP_INCLUDED
#define MATRIX_HPP_INCLUDED

#include "vector.h"

#pragma once
#if !defined(__WIN32) && !defined(__declspec)
#define __declspec(x)
#endif
//Defines the + and - operators for a matrix.
#define DECLARE_FREE_OPERATORS(type) \
	type operator+(const type& rhs, const type& lhs);	\
	type operator-(const type& rhs, const type& lhs);	\
	type operator*(const type& rhs, const type& lhs);	\
	type operator+(float, const type& lhs);				\
	type operator*(float, const type& lhs);				\
	type operator+(const type& lhs, float);				\
	type operator*(const type& lhs, float);				\
	type operator-(const type& lhs, float)

/*
	The following classes implement 2x2, 3x3, and 4x4 square matrices.
	These matrices are stored in row-major order; and can have that data 
	pulled from them by using .data().

	Vectors are treated as column vectors.

	The default constructor constructs the identity matrix.
	
	When indexing a matrix; the upper left corner is 0,0.
	For a 4x4 Matrix; that means that the indexies are

		00 10 20 30
		01 11 21 31
		02 12 22 32
		03 13 23 33

	The following operators are defined.
		+ +=	for both floats and matrices. Component-wise addition.
		-		for floats on the right hand side; and matrices. Component-wise subtraction
		-=		for both floats and matrices. Component-wise subtraction.
		* *=	for floats, matrices and vectors. Note that * for vectors requires the vector
				to be on the right hand side. Component-wise multiplication for floats.
				Matrix multiplication for matrices and vectors.
*/
class Matrix2
{
	friend std::ostream& operator<<(std::ostream&, const Matrix2&);
public:
	static Matrix2 identity();
	static Matrix2 rotation(float amount);
	static Matrix2 scale(const Vector2& amount);
	static Matrix2 scale(float x, float y);
	static Matrix2 scale(float amount);
	static Matrix2 reflection(const Vector2& axis);

	/*
	This function creates the matrix
	x1 x2
	y1 y2

	*/
	Matrix2(float x1, float x2, float y1, float y2);
	Matrix2(const Vector2& col1, const Vector2& col2);
	Matrix2();

	float& operator()(size_t x, size_t y);
	const float& operator()(size_t x, size_t y) const;

	//Component math
	Matrix2& operator+=(const Matrix2&);
	Matrix2& operator-=(const Matrix2&);
	Matrix2& operator+=(float);
	Matrix2& operator-=(float);
	Matrix2	 componentMultiply(const Matrix2&) const;
	Matrix2& componentMultiplyAndAssign(const Matrix2&);
	Matrix2  componentMultiply(float) const;
	Matrix2& componentMultiplyAndAssign(float);

	//A scale.
	Matrix2& operator*=(float);

	//Matrix multiplication
	Matrix2& operator*=(const Matrix2&);

	//Comparison
	bool	operator==(const Matrix2&) const;
	bool	operator!=(const Matrix2&) const;
	
	//Other Math
	float getDeterminant()  const;
	void invert();
	const float * data() const;
private:
	void set(const float x1, const float y1, const float x2, const float y2);
	void set(const Vector2& col1, const Vector2& col2);

	float components[4];
};

DECLARE_FREE_OPERATORS(Matrix2);
Vector2 operator*(const Vector2& row, const Matrix2& m);

class Matrix3
{
	friend  std::ostream& operator<<(std::ostream&, const Matrix3&);
public:
	static Matrix3 rotation(const Vector3& radians);
	static Matrix3 rotation(float x, float y, float z);
	static Matrix3 scale(const Vector3& amount);
	static Matrix3 scale(float x, float y, float z);
	static Matrix3 scale(float amount);
	static Matrix3 identity();
	static Matrix3 face(const Vector3& dir, const Vector3& up);

	//Reflection against the plane m[0]*x + m[1]*y + m[2]*z = 0
	static Matrix3 reflection( const Vector3& m );

	Matrix3(float x1, float x2, float x3,
			float y1, float y2, float y3,
			float z1, float z2, float z3);

	Matrix3(const Vector3& col1, const Vector3& col2, const Vector3& col3);
	Matrix3();

	
	float& operator()(size_t x, size_t y);
	const float& operator()(size_t x, size_t y) const;

	//Component math
	Matrix3& operator+=(const Matrix3&);
	Matrix3& operator-=(const Matrix3&);
	Matrix3& operator+=(const float);
	Matrix3& operator-=(const float);
	Matrix3	 componentMultiply(const Matrix3&) const;
	Matrix3& componentMultiplyAndAssign(const Matrix3&);
	Matrix3  componentMultiply(float) const;
	Matrix3& componentMultiplyAndAssign(float);

	//Multiplication
	Matrix3& operator*=(const float);
	Matrix3& operator*=(const Matrix3&);

	//Comparison
	bool	operator==(const Matrix3&) const;
	bool	operator!=(const Matrix3&) const;

	//Other Math
	float	getDeterminant()  const;
	void invert();
	const float * data() const;
private:
	void set(float x1, float x2, float x3,
		float y1, float y2, float y3,
		float z1, float z2, float z3);
	void set(const Vector3& col1, const Vector3& col2, const Vector3& col3);

	float components[9];
};

DECLARE_FREE_OPERATORS(Matrix3);
Vector3 operator*(const Vector3& row, const Matrix3& m);

//This is the only matrix that you really need.
class Matrix4
{
public:
	static Matrix4 rotation(const Vector4& radians);
	static Matrix4 rotation(float x, float y, float z);
	static Matrix4 scale(const Vector4& amount);
	static Matrix4 scale(float x, float y, float z);
	static Matrix4 scale(float amount);
	static Matrix4 identity();
	static Matrix4 translate(float x, float y, float z);
	static Matrix4 translate(const Vector4& amount);
	static Matrix4 reflection(const Vector4& m);
	static Matrix4 face(const Vector4& dir, const Vector4& up);
	
	friend  std::ostream& operator<<(std::ostream&, const Matrix4&);

	Matrix4(float x1, float x2, float x3, float x4,
			float y1, float y2, float y3, float y4,
			float z1, float z2, float z3, float z4,
			float w1, float w2, float w3, float w4);
	
	Matrix4(const Vector4& col1, const Vector4& col2, const Vector4& col3, const Vector4& col4);
	Matrix4();

	float& operator()(size_t x, size_t y);
	const float& operator()(size_t x, size_t y) const;

	//Component math
	Matrix4& operator+=(const Matrix4&);
	Matrix4& operator-=(const Matrix4&);
	Matrix4& operator+=(float);
	Matrix4& operator-=(float);
	Matrix4	 componentMultiply(const Matrix4&) const;
	Matrix4& componentMultiplyAndAssign(const Matrix4&);
	Matrix4	 componentMultiply(float) const;
	Matrix4& componentMultiplyAndAssign(float);

	//Multiplication
	Matrix4& operator*=(float);
	Matrix4& operator*=(const Matrix4&);
	
	//Comparison
	bool	operator==(const Matrix4&) const;
	bool	operator!=(const Matrix4&) const;

	//Other Math
	float getDeterminant() const;
	void invert();
	const float * data() const;
private:
	void set(float x1, float x2, float x3, float x4,
		float y1, float y2, float y3, float y4,
		float z1, float z2, float z3, float z4,
		float w1, float w2, float w3, float w4);
	void set(const Vector4& row1, const Vector4& row2, const Vector4& row3, const Vector4& row4);
	float components[16];
};

DECLARE_FREE_OPERATORS(Matrix4);
Vector4 operator*(const Vector4& row, const Matrix4& m);

//A template to allow for compile-time matrix size deduction.
//Anything that is not a valid size will result in a compile-time error.
template<size_t N>
struct Matrix
{};

template<>
struct Matrix<2>
{
	typedef Matrix2 type;
};

template<>
struct Matrix<3>
{
	typedef Matrix3 type;
};

template<>
struct Matrix<4>
{
	typedef Matrix4 type;
};

Matrix3 augumentMatrix3(const Matrix2&);
Matrix4 augumentMatrix4(const Matrix2&);
Matrix4 augumentMatrix4(const Matrix3&);

template<size_t N>
Math::Vector<float, N> transform(const Math::Vector<float, N>& v, const typename Matrix<N>::type& t)
{
	return v * t;
}

Matrix2 invert(const Matrix2& m);
Matrix3 invert(const Matrix3& m);
Matrix4 invert(const Matrix4& m);

//A quick utility function for setting the translation of a Matrix4
void setTranslation(Matrix4& m, const Vector4& trans);


//Typedefs for aligned math.
typedef __declspec(align(16)) Matrix2 Matrix2A;
typedef __declspec(align(16)) Matrix3 Matrix3A;
typedef __declspec(align(16)) Matrix4 Matrix4A;
#endif
