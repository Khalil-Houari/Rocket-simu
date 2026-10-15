#pragma once

#include "math_tools.hpp"

class quaternion {
private:
    float q0;
    float q1;
    float q2;
    float q3;

public:
    quaternion();
    quaternion(float q0, float q1, float q2, float q3);

    float get_q0() const;
    float get_q1() const;
    float get_q2() const;
    float get_q3() const;

    void normalize();
    M3 rotation() const;
    V3 euler();
    quaternion deriv(V3 angular_vel) const;

    quaternion operator*(float scalar) const;
    quaternion operator+(const quaternion& other) const;
    quaternion operator/(float scalar) const;
};
