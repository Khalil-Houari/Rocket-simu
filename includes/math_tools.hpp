#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>

extern bool state;

// VECOTR 3D


class V3 {
public:
    float x;
    float y;
    float z;

    V3() : x(0.0f), y(0.0f), z(0.0f) {}
    V3(float x, float y, float z) : x(x), y(y), z(z) {}

    V3 operator+(const V3& other) const {
        return V3(x + other.x, y + other.y, z + other.z);
    }
    V3 operator-(const V3& other) const {
    return V3(x - other.x, y - other.y, z - other.z);
    }
    V3 operator*(float scalar) const {
        return V3(x * scalar, y * scalar, z * scalar);
    }
    V3 operator/(float scalar) const {
        return V3(x / scalar, y / scalar, z / scalar);
    }

    float dot(const V3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    float magnitude() const {
        return sqrt(x * x + y * y + z * z);
    }
    V3 crossV_V(const V3& other) const {
        return V3(  other.z*y - other.y*z ,
                    other.x*z - other.z*x ,
                    other.y*x - other.x*y);
    }
};


// MATRIX 3X3

class M3 {
public:
    float m[3][3];

    M3() { std::fill(&m[0][0], &m[0][0] + 9, 0.0f); }

    M3 operator*(const M3& other) const {
        M3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    result.m[i][j] += m[i][k] * other.m[k][j];
                }
            }
        }
        return result;
    }
    M3 operator*(float a ) const {
        M3 result;
        for(int i = 0; i < 3; ++i) {
            for(int j = 0; j<3; ++j) {
                result.m[i][j] = m[i][j] * a;
            }
        }
        return result;
    }
    M3 operator-(const M3& other) const {
        M3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m[i][j] = m[i][j] - other.m[i][j];
            }
        }
        return result;
    }
    M3 operator+(const M3& other) const {
        M3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m[i][j] = m[i][j] + other.m[i][j];
            }
        }
        return result;
    }
    V3 operator*(const V3& vec) const {
        return V3( m[0][0] * vec.x + m[0][1] * vec.y + m[0][2] * vec.z,
                    m[1][0] * vec.x + m[1][1] * vec.y + m[1][2] * vec.z,
                    m[2][0] * vec.x + m[2][1] * vec.y + m[2][2] * vec.z);
    }
    float determinant() const {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
               m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
               m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }
    M3 transpose() const {
        M3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m[i][j] = m[j][i];
            }
        }
        return result;
    }
    M3 coMatrice() const {
        M3 c;
        c.m[0][0] =  (m[1][1] * m[2][2] - m[1][2] * m[2][1]);
        c.m[0][1] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]);
        c.m[0][2] =  (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        c.m[1][0] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]);
        c.m[1][1] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0]);
        c.m[1][2] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]);
        c.m[2][0] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
        c.m[2][1] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]);
        c.m[2][2] =  (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
        return c;
    }
    M3 inverse() const {
        M3 adj ;
        float det = determinant();
        if(std::abs(det) < 1e-6f) {
            std::cout<< "null determinant !";
            state = true;
        }
        float invdet = 1/det;
        adj = coMatrice();
        adj = adj.transpose();
        return (adj*invdet);
    }
};
