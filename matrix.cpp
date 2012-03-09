#include "matrix.h"
#include "vector.h"

#define DEFINE_SINGLE_OPERATOR(ltype, rtype, ret, op) \
	ret operator op(const ltype& lhs, const rtype& rhs) \
	{													\
		ret tmp(lhs);									\
		tmp BOOST_PP_CAT(op, =) rhs;					\
		return tmp;										\
	}

#define DEFINE_SINGLE_OPERATOR_REV(ltype, rtype, ret, op) \
	ret operator op(const ltype& lhs, const rtype& rhs) \
	{													\
		ret tmp(rhs);									\
		tmp BOOST_PP_CAT(op, =) lhs;					\
		return tmp;										\
	}

#define DEFINE_FREE_OPERATORS(type) \
	DEFINE_SINGLE_OPERATOR(type, type, type, +) \
	DEFINE_SINGLE_OPERATOR(type, type, type, -) \
	DEFINE_SINGLE_OPERATOR(type, type, type, *) \
	DEFINE_SINGLE_OPERATOR_REV(float, type, type, +) \
	DEFINE_SINGLE_OPERATOR_REV(float, type, type, *) \
	DEFINE_SINGLE_OPERATOR(type, float, type, +) \
	DEFINE_SINGLE_OPERATOR(type, float, type, *) \
	DEFINE_SINGLE_OPERATOR(type, float, type, -)

//The Matrix2 member functions
#pragma region Matrix2

Matrix2 Matrix2::identity() 
{
	return Matrix2(
		1, 0,
		0, 1);
}

//Requires a unit vector.
/*
	Constructs a matrix in the form: 

	| x^2 - y^2 | 2 * x * y |
	| 2 * x * y | y^2 - x^2 |

	A reflection matrix around (x, y)
	*/
Matrix2 Matrix2::reflection(const Vector2& axis)
{
	return Matrix2(
		axis[0] * axis[0] - axis[1] * axis[1], 2 * axis[0] * axis[1], 
		2 * axis[0] * axis[1], axis[1] * axis[1] - axis[0] * axis[0]);
}

Matrix2 Matrix2::rotation(float radialAmount)
{
	return Matrix2(
		cos(radialAmount), -sin(radialAmount),
		sin(radialAmount), cos(radialAmount));
}

Matrix2 Matrix2::scale(float amount)
{
	return Matrix2(
		amount, 0,
		0, amount);
}

Matrix2 Matrix2::scale(const Vector2& amount)
{
	return Matrix2::scale(amount[0], amount[1]);
}

Matrix2 Matrix2::scale(float x, float y)
{
	return Matrix2(
		x, 0,
		0, y);
}


Matrix2::Matrix2(const float x1, const float y1, const float x2, const float y2)
{
	set(x1, y1,
		x2, y2);
}

Matrix2::Matrix2()
{
	//Identity;
	set(1, 0,
		0, 1);
}

Matrix2::Matrix2(const Vector2& a, const Vector2& b)
{
	set(a, b);
}

void Matrix2::set(const float m11, const float m12, const float m21, const float m22)
{
	components[0] = m11;
	components[1] = m12;
	components[2] = m21;
	components[3] = m22;
}

void Matrix2::set(const Vector2& column1, const Vector2& column2)
{
	set(column1[0], column2[0],
		column1[1], column2[1]);
}

float& Matrix2::operator()(size_t y, size_t x)
{
	assert(x < 2 && y < 2);
	return components[y * 2 + x];
}

const float& Matrix2::operator()(size_t y, size_t x) const
{
	assert(x < 2 && y < 2);
	return components[y * 2 + x];
}

Matrix2& Matrix2::operator+=(const Matrix2& rhv)
{
	components[0] += rhv.components[0];
	components[1] += rhv.components[1];
	components[2] += rhv.components[2];
	components[3] += rhv.components[3];
	return *this;
}

Matrix2& Matrix2::operator-=(const Matrix2& rhv)
{
	components[0] -= rhv.components[0];
	components[1] -= rhv.components[1];
	components[2] -= rhv.components[2];
	components[3] -= rhv.components[3];
	return *this;
}

Matrix2& Matrix2::operator+=(float rhv)
{
	components[0] += rhv;
	components[1] += rhv;
	components[2] += rhv;
	components[3] += rhv;
	return *this;
}

Matrix2& Matrix2::operator-=(float rhv)
{
	components[0] -= rhv;
	components[1] -= rhv;
	components[2] -= rhv;
	components[3] -= rhv;
	return *this;
}

Matrix2& Matrix2::componentMultiplyAndAssign(const Matrix2& rhv) 
{
	components[0] *= rhv.components[0];
	components[1] *= rhv.components[1];
	components[2] *= rhv.components[2];
	components[3] *= rhv.components[3];
	return *this;
}

Matrix2& Matrix2::componentMultiplyAndAssign(float rhv) 
{
	components[0] *= rhv;
	components[1] *= rhv;
	components[2] *= rhv;
	components[3] *= rhv;
	return *this;
}

Matrix2 Matrix2::componentMultiply(const Matrix2& rhv) const 
{
	Matrix2 temp(*this);
	temp.componentMultiplyAndAssign(rhv);
	return temp;
}

Matrix2 Matrix2::componentMultiply(const float rhv) const 
{
	Matrix2 temp(*this);
	temp.componentMultiplyAndAssign(rhv);
	return temp;
}

Matrix2& Matrix2::operator*=(float rhv) 
{
	return componentMultiplyAndAssign(rhv);
}

bool Matrix2::operator==(const Matrix2& rhv) const
{
	float totalEpsilon = abs((components[0] - rhv.components[0])) +
						 abs((components[1] - rhv.components[1])) + 
						 abs((components[2] - rhv.components[2])) +
						 abs((components[3] - rhv.components[3]));
	return (totalEpsilon < 0.0000001);
}

bool Matrix2::operator!=(const Matrix2& rhv) const 
{
	return !((*this) == rhv);
}

Matrix2& Matrix2::operator*=(const Matrix2& rhv)
{
	float m11 = components[0] * rhv.components[0] + components[1] * rhv.components[2];
	float m21 = components[0] * rhv.components[1] + components[1] * rhv.components[3];

	float m12 = components[2] * rhv.components[0] + components[3] * rhv.components[2];
	float m22 = components[2] * rhv.components[1] + components[3] * rhv.components[3];

	set(m11, m21, 
		m12, m22);

	return *this;
}

float Matrix2::getDeterminant() const
{
	return components[0] * components[3] - components[1] * components[2];
}

void Matrix2::invert()
{
	/*
		2x2 Matrix inversion is done with the formula
		
		[ a b ]
		[ c d ]
		
		1 / det * [d -b][-c a]
	*/

	assert(fabs(getDeterminant()) > 0.000001f);
	float invDet = 1 / getDeterminant();

	float a = components[0];
	float b = components[1];
	float c = components[2];
	float d = components[3];

	set( d, -b,
		-c,  a);

	componentMultiplyAndAssign(invDet);
}

Vector2 operator *(const Vector2& rhv, const Matrix2& m)
{
	float x = rhv[0] * m(0, 0) + rhv[1] * m(1, 0);
	float y = rhv[0] * m(0, 1) + rhv[1] * m(1, 1);

	return Vector2(x, y);
}

DEFINE_FREE_OPERATORS(Matrix2)
#pragma endregion

//Matrix3 member functions
#pragma region Matrix3

Matrix3::Matrix3(float x1, float x2, float x3,
				float y1, float y2, float y3,
				float z1, float z2, float z3)
{
	set(x1, x2, x3,
		y1, y2, y3,
		z1, z2, z3);
}

Matrix3::Matrix3()
{
	set(1, 0, 0,
		0, 1, 0,
		0, 0, 1);
}

Matrix3::Matrix3(const Vector3& a, const Vector3& b, const Vector3& c)
{
	set(a, b, c);
}

void Matrix3::set(float x1, float x2, float x3,
				float y1, float y2, float y3,
				float z1, float z2, float z3)
{
	components[0] = x1;
	components[1] = x2;
	components[2] = x3;
			   
	components[3] = y1;
	components[4] = y2;
	components[5] = y3;
			   
	components[6] = z1;
	components[7] = z2;
	components[8] = z3;
}

void Matrix3::set(const Vector3& col1, const Vector3& col2, const Vector3& col3)
{
	set(col1[0], col2[0], col3[0], 
		col1[1], col2[1], col3[1], 
		col1[2], col2[2], col3[2]);
}

float& Matrix3::operator ()(size_t y, size_t x)
{
	assert(x < 3 && y < 3);
	return components[y * 3 + x];
}

const float& Matrix3::operator ()(size_t y, size_t x) const
{
	assert(x < 3 && y < 3);
	return components[y * 3 + x];
}

Matrix3 Matrix3::componentMultiply(const Matrix3& rhv) const
{
	Matrix3 temp(*this);
	temp.componentMultiplyAndAssign(rhv);
	return temp;
}

Matrix3 Matrix3::componentMultiply(float rhv) const
{
	Matrix3 temp(*this);
	temp.componentMultiplyAndAssign(rhv);
	return temp;
}

Matrix3& Matrix3::operator+=(const Matrix3& rhv)
{
	components[0] += rhv.components[0];
	components[1] += rhv.components[1];
	components[2] += rhv.components[2];
						
	components[3] += rhv.components[3];
	components[4] += rhv.components[4];
	components[5] += rhv.components[5];
						
	components[6] += rhv.components[6];
	components[7] += rhv.components[7];
	components[8] += rhv.components[8];

	return *this;
}

Matrix3& Matrix3::operator-=(const Matrix3& rhv)
{
	components[0] -= rhv.components[0];
	components[1] -= rhv.components[1];
	components[2] -= rhv.components[2];
				  
	components[3] -= rhv.components[3];
	components[4] -= rhv.components[4];
	components[5] -= rhv.components[5];
				  
	components[6] -= rhv.components[6];
	components[7] -= rhv.components[7];
	components[8] -= rhv.components[8];

	return *this;
}

Matrix3& Matrix3::componentMultiplyAndAssign(const Matrix3& rhv)
{
	components[0] *= rhv.components[0];
	components[1] *= rhv.components[1];
	components[2] *= rhv.components[2];
				  
	components[3] *= rhv.components[3];
	components[4] *= rhv.components[4];
	components[5] *= rhv.components[5];
				  
	components[6] *= rhv.components[6];
	components[7] *= rhv.components[7];
	components[8] *= rhv.components[8];
	return *this;
}

Matrix3& Matrix3::operator+=(const float rhv)
{
	components[0] += rhv;
	components[1] += rhv;
	components[2] += rhv;

	components[3] += rhv;
	components[4] += rhv;
	components[5] += rhv;

	components[6] += rhv;
	components[7] += rhv;
	components[8] += rhv;

	return *this;
}

Matrix3& Matrix3::operator-=(const float rhv)
{
	components[0] -= rhv;
	components[1] -= rhv;
	components[2] -= rhv;

	components[3] -= rhv;
	components[4] -= rhv;
	components[5] -= rhv;

	components[6] -= rhv;
	components[7] -= rhv;
	components[8] -= rhv;

	return *this;
}

Matrix3& Matrix3::componentMultiplyAndAssign(const float rhv)
{
	components[0] *= rhv;
	components[1] *= rhv;
	components[2] *= rhv;

	components[3] *= rhv;
	components[4] *= rhv;
	components[5] *= rhv;

	components[6] *= rhv;
	components[7] *= rhv;
	components[8] *= rhv;

	return (*this);
}

Matrix3& Matrix3::operator*=(const float rhv)
{
	return componentMultiplyAndAssign(rhv);
}

Matrix3& Matrix3::operator*=(const Matrix3 &rhv)
{
	const float * c = components;
	const float * o = rhv.components;

	float x1 = c[0] * o[0] + c[1] * o[3] + c[2] * o[6];
	float x2 = c[0] * o[1] + c[1] * o[4] + c[2] * o[7];
	float x3 = c[0] * o[2] + c[1] * o[5] + c[2] * o[8];

	float y1 = c[3] * o[0] + c[4] * o[3] + c[5] * o[6];
	float y2 = c[3] * o[1] + c[4] * o[4] + c[5] * o[7];
	float y3 = c[3] * o[2] + c[4] * o[5] + c[5] * o[8];

	float z1 = c[6] * o[0] + c[7] * o[3] + c[8] * o[6];
	float z2 = c[6] * o[1] + c[7] * o[4] + c[8] * o[7];
	float z3 = c[6] * o[2] + c[7] * o[5] + c[8] * o[8];
	
	set(x1, x2, x3, y1, y2, y3, z1, z2, z3);
	return (*this);
}

bool Matrix3::operator==(const Matrix3& rhv) const
{
	float totalEpsilon = 
		abs((components[0] - rhv.components[0])) + 
		abs((components[1] - rhv.components[1])) + 
		abs((components[2] - rhv.components[2])) + 
		abs((components[3] - rhv.components[3])) + 
		abs((components[4] - rhv.components[4])) +
		abs((components[5] - rhv.components[5])) +
		abs((components[6] - rhv.components[6])) +
		abs((components[7] - rhv.components[7])) +
		abs((components[8] - rhv.components[8]));
	return (totalEpsilon < 0.00000001);
}

bool Matrix3::operator!=(const Matrix3& rhv) const
{
	return !((*this) == rhv);
}

float Matrix3::getDeterminant() const
{
	const float * c = components;
	float det = 
		c[0] * c[4] * c[8] + c[1] * c[5] * c[6] + c[2] * c[3] * c[7] - 
		c[2] * c[4] * c[6] - c[1] * c[3] * c[8] - c[0] * c[5] * c[7];
	return det;
}

void Matrix3::invert()
{
	//See http://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3.C3.973_matrices
	//for the algorithm used here to invert 3x3 matrices.
	float det = getDeterminant();
	assert(fabs(det) > 0.000001f);

	float invdet = 1 / det;

	//With
	//a=0 b=1 c=2
	//d=3 e=4 f=5
	//g=6 h=7 k=8

	//ek - fh
	float a = components[4] * components[8] - components[5] * components[7];

	//ch - bk
	float b = components[2] * components[7] - components[1] * components[8];
	
	//bf - ce
	float c = components[1] * components[5] - components[2] * components[4];

	//fg - dk
	float d = components[5] * components[6] - components[3] * components[8];

	//ak - cg
	float e = components[0] * components[8] - components[2] * components[6];

	//cd - af
	float f = components[2] * components[3] - components[0] * components[5];

	//dh - eg
	float g = components[3] * components[7] - components[4] * components[6];

	//bg - ah
	float h = components[1] * components[6] - components[0] * components[7];

	//ae - bd
	float k = components[0] * components[4] - components[1] * components[3];

	set(a, b, c, 
		d, e, f, 
		g, h, k);

	componentMultiplyAndAssign(invdet);
}

const float * Matrix3::data() const
{
	return components;
}

Matrix3 Matrix3::identity()
{
	Matrix3 id(1, 0, 0,
			   0, 1, 0,
			   0, 0, 1);
	return id;
}

Matrix3 Matrix3::reflection(const Vector3& in)
{
	float a = in[0];
	float b = in[1];
	float c = in[2];
	return Matrix3
		(1 - 2 * a * 1, -2 * a * b, -2 * a * c, 
		 -2 * a * b, 1 - 2 * b * b, -2 * b * c,  
		 -2 * a * c, -2 * b * c, 1 - 2 * c * c);
}

Matrix3 Matrix3::rotation(float x, float y, float z)
{
	Matrix3 xRotation(
		1, 0, 0,
		0, cos(x), sin(x),
		0, -sin(x), cos(x));

	Matrix3 yRotation(
		cos(y), 0, -sin(y),
		0, 1, 0,
		sin(y), 0, cos(y));

	Matrix3 zRotation(
		cos(z), sin(z), 0,
		-sin(z), cos(z), 0,
		0, 0, 1);

	return xRotation * yRotation * zRotation;
}


Matrix3 Matrix3::rotation(const Vector3& radians)
{
	return Matrix3::rotation(radians[0], radians[1], radians[2]);
}

Matrix3 Matrix3::scale(float x, float y, float z)
{
	return Matrix3(
		x, 0, 0,
		0, y, 0,
		0, 0, z);
}

Matrix3 Matrix3::scale(const Vector3& amount)
{
	return Matrix3::scale(amount[0], amount[1], amount[2]);
}

Matrix3 Matrix3::scale(float amount)
{
	return Matrix3::scale(amount, amount, amount);
}

Matrix3 Matrix3::face(const Vector3& d, const Vector3& up)
{
	//Generate a set of new basis vectors and use these.

	Vector3 r = d.perpendicularTo();
	Vector3 u = cross(d, r);

	return Matrix3(d, r, u);
}

Vector3 operator *(const Vector3& rhv, const Matrix3& m)
{
	float x = rhv[0] * m(0, 0) + rhv[1] * m(1, 0) + rhv[2] * m(2, 0);
	float y = rhv[0] * m(0, 1) + rhv[1] * m(1, 1) + rhv[2] * m(2, 1);
	float z = rhv[0] * m(0, 2) + rhv[1] * m(1, 2) + rhv[2] * m(2, 2);

	return Vector3(x, y, z);
}
DEFINE_FREE_OPERATORS(Matrix3)
#pragma endregion


//Matrix 4 member functions.
#pragma region Matrix4

Matrix4 Matrix4::translate(float x, float y, float z)
{
	return Matrix4(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		x, y, z, 1);
}
Matrix4 Matrix4::translate(const Vector4& amount)
{
	return Matrix4::translate(amount[0], amount[1], amount[2]);
}

Matrix4 Matrix4::reflection( const Vector4& in)
{
	float a = in[0];
	float b = in[1];
	float c = in[2];
	return Matrix4
		(1 - 2 * a * 1, -2 * a * b, -2 * a * c, 0,
		 -2 * a * b, 1 - 2 * b * b, -2 * b * c, 0,
		 -2 * a * c, -2 * b * c, 1 - 2 * c * c, 0,
		 0, 0, 0, 1);
}

Matrix4 Matrix4::rotation(float x, float y, float z)
{
	Matrix4 xRotation(
		1, 0, 0, 0, 
		0, cos(x), sin(x), 0,
		0, -sin(x), cos(x), 0, 
		0, 0, 0, 1);

	Matrix4 yRotation(
		cos(y), 0, -sin(y), 0,
		0, 1, 0, 0,
		sin(y), 0, cos(y), 0,
		0, 0, 0, 1);

	Matrix4 zRotation(
		cos(z), sin(z), 0, 0,
		-sin(z), cos(z), 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1);

	return xRotation * yRotation * zRotation;
}

Matrix4 Matrix4::rotation(const Vector4& radians)
{
	return Matrix4::rotation(radians[0], radians[1], radians[2]);
}

Matrix4 Matrix4::scale(float x, float y, float z)
{
	return Matrix4(
		x, 0, 0, 0,
		0, y, 0, 0,
		0, 0, z, 0,
		0, 0, 0, 1);
}

Matrix4 Matrix4::face(const Vector4& d, const Vector4& up)
{
	Vector4 r = d.perpendicularTo();
	Vector4 u = cross(d, r);

	return Matrix4(d, r, u, Vector4(0, 0, 0, 1));
}

Matrix4 Matrix4::scale(const Vector4& amount)
{
	return Matrix4::scale(amount[0], amount[1], amount[2]);
}

Matrix4 Matrix4::scale(float amount)
{
	return Matrix4::scale(amount, amount, amount);
}

Matrix4 Matrix4::identity()
{
	return Matrix4();
}

Matrix4::Matrix4(const float x1, const float x2, const float x3, const float x4
				,const float y1, const float y2, const float y3, const float y4
				,const float z1, const float z2, const float z3, const float z4
				,const float w1, const float w2, const float w3, const float w4)
{
	set(x1, x2, x3, x4, 
		y1, y2, y3, y4, 
		z1, z2, z3, z4, 
		w1, w2, w3, w4);
}

Matrix4::Matrix4(const Vector4& col1, const Vector4& col2, const Vector4& col3, const Vector4 &col4)
{
	set(col1, col2, col3, col4);
}

void Matrix4::set(const float x1, const float x2, const float x3, const float x4,
				const float y1, const float y2, const float y3, const float y4,
				const float z1, const float z2, const float z3, const float z4,
				const float w1, const float w2, const float w3, const float w4)
{
	components[0] = x1;
	components[1] = x2;
	components[2] = x3;
	components[3] = x4;

	components[4] = y1;
	components[5] = y2;
	components[6] = y3;
	components[7] = y4;

	components[8] = z1;
	components[9] = z2;
	components[10] = z3;
	components[11] = z4;

	components[12] = w1;
	components[13] = w2;
	components[14] = w3;
	components[15] = w4;
}

void Matrix4::set(const Vector4 &col1, const Vector4 &col2, const Vector4 &col3, const Vector4 &col4)
{
	set( col1[0], col2[0], col3[0], col4[0],
		 col1[1], col2[1], col3[1], col4[1],
		 col1[2], col2[2], col3[2], col4[2],
		 0, 0, 0, 1);
}

Matrix4::Matrix4()
{
	set(1, 0, 0, 0, 
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1);
}

float& Matrix4::operator ()(size_t y, size_t x)
{
	assert(x < 4 && y < 4);
	return components[y * 4 + x];
}

const float& Matrix4::operator ()(size_t y, size_t x) const
{
	assert(x < 4 && y < 4);
	return components[y * 4 + x];
}

Matrix4& Matrix4::operator+=(const Matrix4& rhv) 
{
	for (size_t i = 0; i < 16; ++i)
	{
		components[i] += rhv.components[i];
	}
	return *this;
}

Matrix4& Matrix4::operator-=(const Matrix4& rhv) 
{
	for (size_t i = 0; i < 16; ++i)
	{
		components[i] -= rhv.components[i];
	}
	return *this;
}

Matrix4& Matrix4::operator+=(float rhv) 
{
	for (size_t i = 0; i < 16; ++i)
	{
		components[i] += rhv;
	}
	return *this;
}

Matrix4& Matrix4::operator-=(float rhv) 
{
	for (size_t i = 0; i < 16; ++i)
	{
		components[i] -= rhv;
	}
	return *this;
}

Matrix4& Matrix4::operator*=(float rhv) 
{
	for (size_t i = 0; i < 16; ++i)
	{
		components[i] *= rhv;
	}
	return *this;
}

Matrix4 Matrix4::componentMultiply(const Matrix4& rhv) const
{
	Matrix4 temp(*this);
	temp.componentMultiplyAndAssign(rhv);
	return temp;
}

Matrix4& Matrix4::componentMultiplyAndAssign(const Matrix4& rhv)
{
	for (size_t i = 0; i < 16; ++i)
	{
		components[i] *= rhv.components[i];
	}
	return *this;
}

Matrix4 Matrix4::componentMultiply(float rhv) const
{
	Matrix4 temp(*this);
	temp.componentMultiplyAndAssign(rhv);
	return temp;
}

Matrix4& Matrix4::componentMultiplyAndAssign(float rhv)
{
	for (size_t i = 0; i < 16; ++i)
	{
		components[i] *= rhv;
	}
	return *this;
}

bool Matrix4::operator==(const Matrix4& rhv) const
{
	float totalEpsilon = 
		abs((components[0]  - rhv.components[0] )) + 
		abs((components[1]  - rhv.components[1] )) + 
		abs((components[2]  - rhv.components[2] )) + 
		abs((components[3]  - rhv.components[3] )) + 
		abs((components[4]  - rhv.components[4] )) + 
		abs((components[5]  - rhv.components[5] )) +
		abs((components[6]  - rhv.components[6] )) +
		abs((components[7]  - rhv.components[7] )) +
		abs((components[8]  - rhv.components[8] )) +
		abs((components[9]  - rhv.components[9] )) +
		abs((components[10] - rhv.components[10])) + 
		abs((components[11] - rhv.components[11])) +
		abs((components[12] - rhv.components[12])) +
		abs((components[13] - rhv.components[13])) +
		abs((components[14] - rhv.components[14])) + 
		abs((components[15] - rhv.components[15]));
	return (totalEpsilon < 0.0000001);
}

bool Matrix4::operator!=(const Matrix4& rhv) const
{
	return !((*this)==rhv);
}

Matrix4& Matrix4::operator*=(const Matrix4& rhv)
{
	const float * c = components;
	const float * o = rhv.components;

	float x1 = c[0] * o[0] + c[1] * o[4] + c[2] * o[8] +  c[3] * o[12];
	float x2 = c[0] * o[1] + c[1] * o[5] + c[2] * o[9] +  c[3] * o[13];
	float x3 = c[0] * o[2] + c[1] * o[6] + c[2] * o[10] + c[3] * o[14];
	float x4 = c[0] * o[3] + c[1] * o[7] + c[2] * o[11] + c[3] * o[15];
	
	float y1 = c[4] * o[0] + c[5] * o[4] + c[6] * o[8] +  c[7] * o[12];
	float y2 = c[4] * o[1] + c[5] * o[5] + c[6] * o[9] +  c[7] * o[13];
	float y3 = c[4] * o[2] + c[5] * o[6] + c[6] * o[10] + c[7] * o[14];
	float y4 = c[4] * o[3] + c[5] * o[7] + c[6] * o[11] + c[7] * o[15];

	float z1 = c[8] * o[0] + c[9] * o[4] + c[10] * o[8] +  c[11] * o[12];
	float z2 = c[8] * o[1] + c[9] * o[5] + c[10] * o[9] +  c[11] * o[13];
	float z3 = c[8] * o[2] + c[9] * o[6] + c[10] * o[10] + c[11] * o[14];
	float z4 = c[8] * o[3] + c[9] * o[7] + c[10] * o[11] + c[11] * o[15];

	float w1 = c[12] * o[0] + c[13] * o[4] + c[14] * o[8] +  c[15] * o[12];
	float w2 = c[12] * o[1] + c[13] * o[5] + c[14] * o[9] +  c[15] * o[13];
	float w3 = c[12] * o[2] + c[13] * o[6] + c[14] * o[10] + c[15] * o[14];
	float w4 = c[12] * o[3] + c[13] * o[7] + c[14] * o[11] + c[15] * o[15];

	set(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, w1, w2, w3, w4);
	return (*this);
}

//Code from geometrictools.com
float Matrix4::getDeterminant() const
{
	float a0 = components[ 0] * components[ 5] - components[ 1] * components[ 4];
	float a1 = components[ 0] * components[ 6] - components[ 2] * components[ 4];
	float a2 = components[ 0] * components[ 7] - components[ 3] * components[ 4];
	float a3 = components[ 1] * components[ 6] - components[ 2] * components[ 5];
	float a4 = components[ 1] * components[ 7] - components[ 3] * components[ 5];
	float a5 = components[ 2] * components[ 7] - components[ 3] * components[ 6];
	float b0 = components[ 8] * components[13] - components[ 9] * components[12];
	float b1 = components[ 8] * components[14] - components[10] * components[12];
	float b2 = components[ 8] * components[15] - components[11] * components[12];
	float b3 = components[ 9] * components[14] - components[10] * components[13];
	float b4 = components[ 9] * components[15] - components[11] * components[13];
	float b5 = components[10] * components[15] - components[11] * components[14];

	return a0*b5 - a1*b4 + a2*b3 + a3*b2 - a4*b1 + a5*b0;
}

const float * Matrix4::data() const
{
	return components;
}

void Matrix4::invert()
{
	float a0 = components[ 0] * components[ 5] - components[ 1] * components[ 4];
	float a1 = components[ 0] * components[ 6] - components[ 2] * components[ 4];
	float a2 = components[ 0] * components[ 7] - components[ 3] * components[ 4];
	float a3 = components[ 1] * components[ 6] - components[ 2] * components[ 5];
	float a4 = components[ 1] * components[ 7] - components[ 3] * components[ 5];
	float a5 = components[ 2] * components[ 7] - components[ 3] * components[ 6];
	float b0 = components[ 8] * components[13] - components[ 9] * components[12];
	float b1 = components[ 8] * components[14] - components[10] * components[12];
	float b2 = components[ 8] * components[15] - components[11] * components[12];
	float b3 = components[ 9] * components[14] - components[10] * components[13];
	float b4 = components[ 9] * components[15] - components[11] * components[13];
	float b5 = components[10] * components[15] - components[11] * components[14];

	float det = a0*b5 - a1*b4 + a2*b3 + a3*b2 - a4*b1 + a5*b0;
	assert(fabs(det) > 0.00001f);

	float result[16];
	result[ 0] = + components[ 5]*b5 - components[ 6]*b4 + components[ 7]*b3;
	result[ 1] = - components[ 1]*b5 + components[ 2]*b4 - components[ 3]*b3;
	result[ 2] = + components[13]*a5 - components[14]*a4 + components[15]*a3;
	result[ 3] = - components[ 9]*a5 + components[10]*a4 - components[11]*a3;
	result[ 4] = - components[ 4]*b5 + components[ 6]*b2 - components[ 7]*b1;
	result[ 5] = + components[ 0]*b5 - components[ 2]*b2 + components[ 3]*b1;
	result[ 6] = - components[12]*a5 + components[14]*a2 - components[15]*a1;
	result[ 7] = + components[ 8]*a5 - components[10]*a2 + components[11]*a1;
	result[ 8] = + components[ 4]*b4 - components[ 5]*b2 + components[ 7]*b0;
	result[ 9] = - components[ 0]*b4 + components[ 1]*b2 - components[ 3]*b0;
	result[10] = + components[12]*a4 - components[13]*a2 + components[15]*a0;
	result[11] = - components[ 8]*a4 + components[ 9]*a2 - components[11]*a0;
	result[12] = - components[ 4]*b3 + components[ 5]*b1 - components[ 6]*b0;
	result[13] = + components[ 0]*b3 - components[ 1]*b1 + components[ 2]*b0;
	result[14] = - components[12]*a3 + components[13]*a1 - components[14]*a0;
	result[15] = + components[ 8]*a3 - components[ 9]*a1 + components[10]*a0;


	set(result[0], result[1], result[2], result[3], 
		result[4], result[5], result[6], result[7], 
		result[8], result[9], result[10], result[11], 
		result[12], result[13], result[14], result[15]);

	componentMultiplyAndAssign(1 / det);
}

Vector4 operator*(const Vector4& rhv, const Matrix4& m)
{
	float x = rhv[0] * m(0, 0) + rhv[1] * m(1, 0) + rhv[2] * m(2, 0) + rhv[3] * m(3, 0);
	float y = rhv[0] * m(0, 1) + rhv[1] * m(1, 1) + rhv[2] * m(2, 1) + rhv[3] * m(3, 1);
	float z = rhv[0] * m(0, 2) + rhv[1] * m(1, 2) + rhv[2] * m(2, 2) + rhv[3] * m(3, 2);
	float w = rhv[0] * m(0, 3) + rhv[1] * m(1, 3) + rhv[2] * m(2, 3) + rhv[3] * m(3, 3);

	return Vector4(x, y, z, w);
}

DEFINE_FREE_OPERATORS(Matrix4);
#pragma endregion


#pragma region Free Matrix Functions
std::ostream& operator<<(std::ostream& o, const Matrix2& rhv)
{
	o<<"[["<<rhv(0,0)<<","<<rhv(0,1)<<"]["<<rhv(1,0)<<","<<rhv(1,1)<<"]]";
	return o;
}

std::ostream& operator<<(std::ostream& o, const Matrix3& rhv)
{
	o<<"[["<<rhv(0,0)<<","<<rhv(0,1)<<","<<rhv(0,2)<<"]["<<
			 rhv(1,0)<<","<<rhv(1,1)<<","<<rhv(1,2)<<"]["<<
			 rhv(2,0)<<","<<rhv(2,1)<<","<<rhv(2,2)<<"]]";
	return o;
}

std::ostream& operator<<(std::ostream& o, const Matrix4& rhv)
{

	o<<"[["<<rhv(0,0)<<","<<rhv(0,1)<<","<<rhv(0,2)<<","<<rhv(0,3)<<"]["<<
		     rhv(1,0)<<","<<rhv(1,1)<<","<<rhv(1,2)<<","<<rhv(1,3)<<"]["<<
			 rhv(2,0)<<","<<rhv(2,1)<<","<<rhv(2,2)<<","<<rhv(2,3)<<"]["<<
			 rhv(3,0)<<","<<rhv(3,1)<<","<<rhv(3,2)<<","<<rhv(3,3)<<"]]";
	return o;
}

Matrix3 augumentMatrix3(const Matrix2& rhv)
{
	Matrix3 result(rhv(0, 0), rhv(0, 1), 0,
				   rhv(1, 0), rhv(1, 1), 0,
				   0		, 0		   , 1);
	return result;
}

Matrix4 augumentMatrix4(const Matrix3& rhv)
{
	Matrix4 result(rhv(0, 0), rhv(0, 1), rhv(0, 2), 0, 
				   rhv(1, 0), rhv(1, 1), rhv(1, 2), 0,
				   rhv(2, 0), rhv(2, 1), rhv(2, 2), 0,
				   0		, 0		   , 0		  , 1);
	return result;
}

Matrix4 augumentMatrix4(const Matrix2& rhv)
{
	Matrix4 result(
		rhv(0, 0), rhv(0, 1), 0, 0,
		rhv(1, 0), rhv(1, 1), 0, 0,
		0        , 0        , 1, 0,
		0        , 0        , 0, 1);
	return result;
}

Matrix2 invert(const Matrix2& m)
{
	Matrix2 temp(m);
	temp.invert();
	return temp;
}

Matrix3 invert(const Matrix3& m)
{
	Matrix3 temp(m);
	temp.invert();
	return temp;
}

Matrix4 invert(const Matrix4& m)
{
	Matrix4 temp(m);
	temp.invert();
	return temp;
}


void setTranslation(Matrix4& m, const Vector4& trans)
{
	m(3, 0) = trans[0];
	m(3, 1) = trans[1];
	m(3, 2) = trans[2];
}
#pragma endregion
