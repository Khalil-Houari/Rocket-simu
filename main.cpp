#include <cstdio>
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <algorithm>


float dt; // will be defined , it is not the time step for RK4 

class quaternion {
private:
    float q0;
    float q1;
    float q2;
    float q3;

public:
    quaternion() : q0(1.0f), q1(0.0f), q2(0.0f), q3(0.0f) {}

    quaternion(float q0, float q1, float q2, float q3) : q0(q0), q1(q1), q2(q2), q3(q3) {}

    float get_q0() const { return q0; }
    float get_q1() const { return q1; }
    float get_q2() const { return q2; }
    float get_q3() const { return q3; }

    void normalize() {
        float norm = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
        q0 /= norm;
        q1 /= norm;
        q2 /= norm;
        q3 /= norm;
    }

    void rotation(float R[3][3]) const { // rotation matrices
        //i 
        R[0][0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
        R[0][1] = 2 * (q1 * q2 - q0 * q3);
        R[0][2] = 2 * (q1 * q3 + q0 * q2);
        //j
        R[1][0] = 2 * (q1 * q2 + q0 * q3);
        R[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
        R[1][2] = 2 * (q2 * q3 - q0 * q1);
        //k
        R[2][0] = 2 * (q1 * q3 - q0 * q2);
        R[2][1] = 2 * (q2 * q3 + q0 * q1);
        R[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
    }

    void euler(float& alpha, float& beta, float& gamma) {
        alpha = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2));
        
        float sin_beta = 2 * (q0 * q2 - q3 * q1);
        beta = asin(std::clamp(sin_beta, -0.999999f, 0.999999f)); 

        gamma = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3));
    }

    quaternion deriv(float wx, float wy, float wz) const {
        float q0Dot = 0.5f * (-q1 * wx - q2 * wy - q3 * wz);
        float q1Dot = 0.5f * (q0 * wx + q2 * wz - q3 * wy);
        float q2Dot = 0.5f * (q0 * wy - q1 * wz + q3 * wx);
        float q3Dot = 0.5f * (q0 * wz + q1 * wy - q2 * wx);
        return { q0Dot, q1Dot, q2Dot, q3Dot };
    }

    quaternion update(float wx, float wy, float wz, float dt) const {
        quaternion qDot = deriv(wx, wy, wz);
        return quaternion{ q0 + qDot.q0 * dt, q1 + qDot.q1 * dt, q2 + qDot.q2 * dt, q3 + qDot.q3 * dt };
    }

    quaternion operator*(float scalar) const {
        return quaternion{ q0 * scalar, q1 * scalar, q2 * scalar, q3 * scalar };
    }
    quaternion operator+(const quaternion& other) const {
        return quaternion{ q0 + other.q0 , q1 + other.q1, q2 + other.q2, q3 + other.q3 };
    }
    quaternion operator/(float scalar) const {
        return quaternion{ q0 / scalar, q1 / scalar, q2 / scalar, q3 / scalar };
    }
};


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
};

class M3 {
    public:
        float m[3][3];

        M3() {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    m[i][j] = 0.0f;
                }
            }
        }

        M3(float m[3][3]) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    this->m[i][j] = m[i][j];
                }
            }
        }
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
};

class MV {
    public: 
    MV(V3 v, V3 u, M3 m, M3 w) : vec(v), mat(m), vec2(u), mat2(w) {}
    V3 vec;
    M3 mat;
    V3 vec2;
    M3 mat2;

    float dotV_V() {
        return vec.x * vec2.x + vec.y * vec2.y + vec.z * vec2.z;
    }
    V3 crossV_V() {
        return V3(vec.y * vec2.z - vec.z * vec2.y,
                  vec.z * vec2.x - vec.x * vec2.z,
                  vec.x * vec2.y - vec.y * vec2.x);
    }
    M3 mat_mat() {
        return mat * mat2;
    }
    V3 mat_vec(int a, int b) {
        // m1*v1 
        V3 result_vec1( mat.m[0][0] * vec.x + mat.m[0][1] * vec.y + mat.m[0][2] * vec.z,
                        mat.m[1][0] * vec.x + mat.m[1][1] * vec.y + mat.m[1][2] * vec.z,
                        mat.m[2][0] * vec.x + mat.m[2][1] * vec.y + mat.m[2][2] * vec.z);
        //m1*v2
        V3 result_vec2( mat.m[0][0] * vec2.x + mat.m[0][1] * vec2.y + mat.m[0][2] * vec2.z,
                        mat.m[1][0] * vec2.x + mat.m[1][1] * vec2.y + mat.m[1][2] * vec2.z,
                        mat.m[2][0] * vec2.x + mat.m[2][1] * vec2.y + mat.m[2][2] * vec2.z);
        // m2*v1
        V3 result_vec3(mat2.m[0][0] * vec.x + mat2.m[0][1] * vec.y + mat2.m[0][2] * vec.z,
                       mat2.m[1][0] * vec.x + mat2.m[1][1] * vec.y + mat2.m[1][2] * vec.z,
                       mat2.m[2][0] * vec.x + mat2.m[2][1] * vec.y + mat2.m[2][2] * vec.z);
        // m2*v2
        V3 result_vec4(mat2.m[0][0] * vec2.x +  mat2.m[0][1] * vec2.y + mat2.m[0][2] * vec2.z,
                       mat2.m[1][0] * vec2.x +  mat2.m[1][1] * vec2.y + mat2.m[1][2] * vec2.z,
                       mat2.m[2][0] * vec2.x + mat2.m[2][1] * vec2.y + mat2.m[2][2] * vec2.z);
        if (a == 1 && b == 1) {
            return V3(result_vec1);
        } else if (a == 1 && b == 2) {
            return V3(result_vec2);
        } else if (a == 2 && b == 1) {
            return V3(result_vec3);
        } else if (a == 2 && b == 2) {
            return V3(result_vec4);
        } else {
            throw std::invalid_argument("Invalid values for a and b, 1 or 2"); //debugging 
        }
    }
};







// the rk4 is a vestige of a past failed attempt we ll reuse it 



template<typename T>
    class RK4 {
    public:
        static void rk4(R_state& s, const R_dynamics& dynamics, float dt) {
            R_state k1 = dynamics.deriv(s);
            R_state k2 = dynamics.deriv(s + k1 * dt / 2.0f);
            R_state k3 = dynamics.deriv(s + k2 * dt / 2.0f);
            R_state k4 = dynamics.deriv(s + k3 * dt);

            s = s + (k1 + k2 * 2.0f + k3 * 2.0f + k4) * dt / 6.0f;

            s.ori.normalize();

            s.update_ang();
        }
    };





    //reworked all the algebra tools and some mathematical tools 
    // need to implement the rest of the parameters and the dynamics of the rocket




