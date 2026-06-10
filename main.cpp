#include <cstdio>
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <algorithm>
#define PI 3.14159265358979323846f

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
    MV() {}

    float dotV_V(V3 vec, V3 vec2) {
        return vec.x * vec2.x + vec.y * vec2.y + vec.z * vec2.z;
    }
    V3 crossV_V(V3 vec, V3 vec2) {
        return V3(vec.y * vec2.z - vec.z * vec2.y,
                  vec.z * vec2.x - vec.x * vec2.z,
                  vec.x * vec2.y - vec.y * vec2.x);
    }
    M3 mat_mat(M3 mat, M3 mat2) {
        return mat * mat2;
    }
    V3 mat_vec(M3 mat, V3 vec) {
        // m1*v1 
        V3 result_vec1( mat.m[0][0] * vec.x + mat.m[0][1] * vec.y + mat.m[0][2] * vec.z,
                        mat.m[1][0] * vec.x + mat.m[1][1] * vec.y + mat.m[1][2] * vec.z,
                        mat.m[2][0] * vec.x + mat.m[2][1] * vec.y + mat.m[2][2] * vec.z);
        return result_vec1;
    }
};

// helper function 
M3 inertia_cyclindre(float mass, float radius, float height) {
    M3 inertia;
    inertia.m[0][0] = mass * (3.0f *radius * radius + height * height) / 12.0f; // Ixx = Iyy = (1/12) * m * (3r^2 + h^2)
    inertia.m[1][1] = inertia.m[0][0]; 
    inertia.m[2][2] =  0.5f * mass * radius * radius; // Izz = (1/2) * m * r^2
    return inertia;
}

M3 inertia_parallelepiped(float mass, float length, float width, float height) {
    M3 inertia;
    inertia.m[0][0] = (1.0f / 12.0f) * mass * (width * width + height * height);
    inertia.m[1][1] = (1.0f / 12.0f) * mass * (length * length + height * height);
    inertia.m[2][2] = (1.0f / 12.0f) * mass * (length * length + width * width);
    return inertia;
}

M3 inertia_rhead(float density, float radius, float height) {
    M3 inertia;
    inertia.m[0][0] = PI*density*(32.0f * radius * radius * radius * radius * height/315.0f + 29*radius * radius * height * height * height/105.0f); 
    inertia.m[1][1] = inertia.m[0][0];
    inertia.m[2][2] = 0.5f * PI * 128.0f * density * radius * radius * radius * radius * height/315.0f; // 
    return inertia;
}



struct rocket_parameters { // still sketching 
    //shell parameters
    float mass_shell;
    float mass_head;
    float density_shell;
    float radius_ext;
    float thickness;
    float height;
    float head_height;
    float radius_int = radius_ext - thickness;
    float pos_shell;
    float pos_body;
    float pos_head;
    float Cog_shell = pos_shell + (height) / 2.0f; // assuming the shape is a perfect cylindre , a thing to change 
    // motor parameters
    float radius_motor;
    float motor_length;
    float mass_motor;
    float dm; // mass flow rate of the fuel, this one too
    float pos_motor;
    float CoG_motor = pos_motor + motor_length / 2.0f; // the cog is in the middle since the fuel + motor is assumed to be homogeneous, but later on when we would rework the motor parameters we will change it
    // FC parameters
    float radius_FC = radius_int; // for now we will assume that the FC has the same radius as the inner radius of the shell, but we can change it later 
    float mass_FC;
    float length_FC;
    float pos_FC;
    float CoG_FC = pos_FC + length_FC / 2.0f; // we will assume that the FC is homogeneous
    // nanosat parameters
    float length_nanosat;
    float width_nanosat;
    float height_nanosat;
    float mass_nanosat;
    float pos_nanosat;
    float CoG_nanosat = pos_nanosat + length_nanosat / 2.0f; // we will assume that the nanosat is homogeneous and that it is placed at the end of the rocket
   
    float Cd; // drag coeff
    float Rho; // air density


};

class rocket_state { // still sketching
    public:
    rocket_parameters p;
    V3 position;
    V3 velocity;
    quaternion orientation;
    V3 angular_velocity; 
    M3 steiner(float mass, float d) {
    M3 result;
    result.m[0][0] = mass * d * d; 
    result.m[1][1] = mass * d * d; 
    result.m[2][2] = 0.0f;         
    return result;
    }
    

    float mass_fuel; // we will update this one during the simulation
    float dm; // mass flow rate of the fuel
    float get_mass_total() const {
        return p.mass_shell + mass_fuel + p.mass_motor + p.mass_FC + p.mass_nanosat;
    }

    float get_CoG() const {
        return (p.mass_shell * p.Cog_shell + p.mass_motor * p.CoG_motor 
               + p.mass_FC * p.CoG_FC + p.mass_nanosat * p.CoG_nanosat 
               + mass_fuel * p.CoG_motor) / get_mass_total();
    }
 // we will assume that the fuel is distributed uniformly along the motor length
    
    rocket_state() : position(), velocity(), orientation(), angular_velocity() {}
    rocket_state(V3 pos, V3 vel, quaternion orient, V3 ang_vel) : position(pos), velocity(vel), orientation(orient), angular_velocity(ang_vel) {  
    }



    M3 get_inertia_tensor(){ 
        float CoG = get_CoG() ;
        return inertia_cyclindre(p.mass_shell - p.mass_head, p.radius_ext, p.height-p.head_height)  +steiner(p.mass_shell-p.mass_head,(p.pos_body-CoG)) 

        + inertia_rhead(p.density_shell, p.radius_ext, p.head_height) +steiner(p.mass_head, p.pos_head- CoG)

        + inertia_cyclindre(p.mass_FC, p.radius_FC*2, p.length_FC) +steiner(p.mass_FC,(p.pos_FC-CoG)) 

        + inertia_parallelepiped(p.mass_nanosat, p.length_nanosat, p.width_nanosat, p.height_nanosat) +steiner(p.mass_nanosat,(p.CoG_nanosat - CoG))
        
        + inertia_cyclindre(p.mass_motor,p.radius_motor,p.motor_length) +steiner(p.mass_motor, p.pos_motor-CoG);
    }
};


class Rk4 {
    public:
    static void rk4(quaternion& q, float wx, float wy, float wz, float dt) {
        quaternion k1 = q.deriv(wx, wy, wz);
        quaternion k2 = (q + k1 * dt / 2.0f).deriv(wx, wy, wz);
        quaternion k3 = (q + k2 * dt / 2.0f).deriv(wx, wy, wz);
        quaternion k4 = (q + k3 * dt).deriv(wx, wy, wz);

        q = q + (k1 + k2 * 2.0f + k3 * 2.0f + k4) * (dt / 6.0f);
        q.normalize();
    }
};

// for now we have only RK4 integrator for the rotations , but well surely add one more integrators for the other components of the rocket's dynamics





    



    // reworked all the algebra tools and some mathematical tools 
    // reworked the quaternion class to be more efficient and more accurate
    // reworeked the rk4 integrator
    // finally finished with inertia tensors and parameters structure even tho we're still missing some parameters 
    // overall the code is getting structred
    // need to start with the dynamics of the engine



