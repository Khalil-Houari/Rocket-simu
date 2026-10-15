#include "quaternion.hpp"
#include <algorithm>
#include <cmath>

quaternion::quaternion() : q0(1.0f), q1(0.0f), q2(0.0f), q3(0.0f) {}

quaternion::quaternion(float q0, float q1, float q2, float q3)
    : q0(q0), q1(q1), q2(q2), q3(q3) {}

float quaternion::get_q0() const { return q0; }
float quaternion::get_q1() const { return q1; }
float quaternion::get_q2() const { return q2; }
float quaternion::get_q3() const { return q3; }

void quaternion::normalize() {
    float norm = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q0 /= norm;
    q1 /= norm;
    q2 /= norm;
    q3 /= norm;
}

M3 quaternion::rotation() const {
    M3 R;
    R.m[0][0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    R.m[0][1] = 2 * (q1 * q2 - q0 * q3);
    R.m[0][2] = 2 * (q1 * q3 + q0 * q2);
    R.m[1][0] = 2 * (q1 * q2 + q0 * q3);
    R.m[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    R.m[1][2] = 2 * (q2 * q3 - q0 * q1);
    R.m[2][0] = 2 * (q1 * q3 - q0 * q2);
    R.m[2][1] = 2 * (q0 * q1 + q2 * q3);
    R.m[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
    return R;
}

V3 quaternion::euler() {
    float x = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2));
    float sin_beta = 2 * (q0 * q2 - q3 * q1);
    float y = asin(std::clamp(sin_beta, -0.999999f, 0.999999f));
    float z = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3));
    return V3(x, y, z);
}

quaternion quaternion::deriv(V3 angular_vel) const {
    float q0Dot = 0.5f * (-q1 * angular_vel.x - q2 * angular_vel.y - q3 * angular_vel.z);
    float q1Dot = 0.5f * (q0 * angular_vel.x + q2 * angular_vel.z - q3 * angular_vel.y);
    float q2Dot = 0.5f * (q0 * angular_vel.y - q1 * angular_vel.z + q3 * angular_vel.x);
    float q3Dot = 0.5f * (q0 * angular_vel.z + q1 * angular_vel.y - q2 * angular_vel.x);
    return {q0Dot, q1Dot, q2Dot, q3Dot};
}

quaternion quaternion::operator*(float scalar) const {
    return quaternion{ q0 * scalar, q1 * scalar, q2 * scalar, q3 * scalar };
}

quaternion quaternion::operator+(const quaternion& other) const {
    return quaternion{ q0 + other.q0 , q1 + other.q1, q2 + other.q2, q3 + other.q3 };
}

quaternion quaternion::operator/(float scalar) const {
    return quaternion{ q0 / scalar, q1 / scalar, q2 / scalar, q3 / scalar };
}
