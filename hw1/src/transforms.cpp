#include "transforms.h"

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

namespace CGL {

Vector2D operator*(const Matrix3x3 &m, const Vector2D &v) {
	Vector3D mv = m * Vector3D(v.x, v.y, 1);
	return Vector2D(mv.x / mv.z, mv.y / mv.z);
}

Matrix3x3 translate(float tx, float ty) {
    Matrix3x3 m;
    m(0,0) = 1; m(0,1) = 0; m(0,2) = tx;
    m(1,0) = 0; m(1,1) = 1; m(1,2) = ty;
    m(2,0) = 0; m(2,1) = 0; m(2,2) = 1;
    return m;
}

Matrix3x3 scale(float sx, float sy) {
    Matrix3x3 m;
    m(0,0) = sx; m(0,1) = 0;  m(0,2) = 0;
    m(1,0) = 0;  m(1,1) = sy; m(1,2) = 0;
    m(2,0) = 0;  m(2,1) = 0;  m(2,2) = 1;
    return m;
}

Matrix3x3 rotate(float deg) {
    float rad = deg * M_PI / 180.0f;
    float cosr = cos(rad);
    float sinr = sin(rad);
    return Matrix3x3(
         cosr, -sinr,  0,
         sinr,  cosr,  0,
         0,    0,    1
    );
}


}
