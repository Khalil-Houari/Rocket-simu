#include <cstdio>
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <algorithm>

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

    void rotation(float R[3][3]) const {
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
        beta = asin(std::clamp(sin_beta, -0.999999f, 0.999999f)); // clamp to avoid numerical issues

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

    V3() : x(0), y(0), z(0) {}

    V3(float x, float y, float z) : x(x), y(y), z(z) {}

    void norm() {
        float norm = sqrt(x * x + y * y + z * z);
        x /= norm;
        y /= norm;
        z /= norm;
    }

    V3 cross(const V3& other) const {
        V3 result;
        result.x = y * other.z - z * other.y;
        result.y = z * other.x - x * other.z;
        result.z = x * other.y - y * other.x;
        return result;
    }

    V3 transform(float R[3][3]) const {
        V3 result;
        result.x = R[0][0] * x + R[0][1] * y + R[0][2] * z;
        result.y = R[1][0] * x + R[1][1] * y + R[1][2] * z;
        result.z = R[2][0] * x + R[2][1] * y + R[2][2] * z;
        return result;
    }

    V3 transform_T(float R[3][3]) const {
        V3 result;
        result.x = R[0][0] * x + R[1][0] * y + R[2][0] * z;
        result.y = R[0][1] * x + R[1][1] * y + R[2][1] * z;
        result.z = R[0][2] * x + R[1][2] * y + R[2][2] * z;
        return result;
    }

    float mag() const {
        return sqrt(x * x + y * y + z * z);
    }

    V3 normalized() const {
        float magn = mag();
        if (magn > 0.0001f)
            return { x / magn, y / magn, z / magn };
        return V3{ 0, 0, 0 };
    }

    V3 operator+(const V3& other) const {
        return { x + other.x, y + other.y, z + other.z };
    }

    V3 operator-(const V3& other) const {
        return { x - other.x, y - other.y, z - other.z };
    }

    V3 operator*(float scalar) const {
        return { x * scalar, y * scalar, z * scalar };
    }

    V3 operator/(float scalar) const {
        return { x / scalar, y / scalar, z / scalar };
    }

};


class R_state {
public:
    V3 pos;
    quaternion ori;
    V3 vel;
    V3 ang_vel;

    float mass;
    float car_mass;

    float pitch = 0;
    float yaw=0;
    float roll = 0;

    R_state() {} 

    void update_ang() {
        ori.euler(pitch, yaw, roll);
    }

    R_state operator+(const R_state& other) const {
        R_state result;
        result.pos = pos + other.pos;
        result.ori = ori + other.ori;
        result.vel = vel + other.vel;
        result.ang_vel = ang_vel + other.ang_vel;
        result.mass = mass + other.mass;
        result.car_mass = car_mass + other.car_mass;
        return result;  
    }
    R_state operator*(float scalar) const {
        R_state result;
        result.pos = pos * scalar;
        result.ori = ori * scalar;
        result.vel = vel * scalar;
        result.ang_vel = ang_vel * scalar;
        result.mass = mass * scalar;
        result.pitch = pitch * scalar;
        result.yaw = yaw * scalar;
        result.roll = roll * scalar;
        result.car_mass = car_mass * scalar;
        return result;
    }

    R_state operator/(float scalar) const {
        R_state result;
        result.pos = pos / scalar;
        result.ori = ori / scalar; 
        result.vel = vel / scalar;
            result.ang_vel = ang_vel / scalar;
            result.mass = mass / scalar;
            result.pitch = pitch / scalar;
            result.yaw = yaw / scalar;
            result.roll = roll / scalar;
            result.car_mass = car_mass / scalar;
            return result;
        }
    };


    class R_params {
    public:
        float m_motor = 20;
        float m_NS= 5;
        float m_FC= 2;
        float m_shell=7;
        float m_car= 10; // carburant
        float dm= -0;


        float k= 2; //drag coeff
        float r = 16; //radius internal
        float R = 17; //radius external

        float d_CoP = -0.7; // CoP from CoG
        float r_thrust = 1; // distance from CoG to thrust
        float d_motor= 0.8; // distance from CoG to motor
        float d_NS= 1; // distance from CoG to NS
        float d_FC = 0.3; // distance from CoG to FC

        float Thrust = 1700; //a constant for now

        R_params() {
            m_motor = 20;
            m_NS = 5;
            m_FC = 2;
            m_shell = 7;
            m_car = 10; // carburant
            dm = -0.5;

            k = 0.002; //drag coeff
            r = 16; //radius internal
            R = 17; //radius external

            d_CoP = -0.7; // CoP from CoG
            r_thrust = 1; // distance from CoG to thrust
            d_motor = 0.8; // distance from CoG to motor
            d_NS = 1; // distance from CoG to NS
            d_FC = 0.3; // distance from CoG to FC

            Thrust = 1700; //a constant for now
        }

        

        float inertia(float f_mass) const {
            float I_shell = m_shell * (R * R / 2 + r * r / 2);
            float I_motor = m_motor * (d_motor * d_motor + r * r / 2);
            float I_NS = m_NS * (d_NS * d_NS + r * r / 2);
            float I_FC = m_FC * (d_FC * d_FC + r * r / 2);
            float I_car = f_mass * (d_motor * d_motor + r * r / 2); // Assuming the fuel is located at the same position as the motor and cylindrical 
            return I_shell + I_motor + I_NS + I_FC + I_car ;
        }
        
        float total_mass() const {
            return m_shell + m_motor + m_NS + m_FC + m_car;
        }
    };


    class R_dynamics {
    private:
        const R_params& params;
    public:
        R_dynamics(const R_params& p) : params(p) {}

        R_state deriv( const R_state& s) const {
            float R[3][3];
            s.ori.rotation(R);

            V3 thrust_rkt = { 0, params.Thrust, 0 }; // assuming thrust is in the y direction of the body frame

            V3 v_rkt = s.vel.transform_T(R);

            V3 f_drag = { 0, 0, 0 }; // assuming drag is in the y direction of the body frame

            float mag_v = v_rkt.mag();
            if (mag_v > 0.001f) {
                float f_drag_mag = -params.k * mag_v * mag_v; // drag force magnitude proportional to the square of velocity
                f_drag = v_rkt.normalized() * f_drag_mag;
            }


            V3 F_T = thrust_rkt + f_drag;

            V3 F = F_T.transform(R);

            float T_mass = s.mass ;

            F.y += -9.81f * T_mass;

            // lever 
            V3 r_thrust = { 0.0f, params.r_thrust, 0.0f };
            V3 r_CoP = { 0.0f, params.d_CoP, 0.0f };

            V3 T_thr = r_thrust.cross(thrust_rkt);
            V3 T_drag = r_CoP.cross(f_drag);
            V3 T = T_thr + T_drag;

            float I = params.inertia(s.car_mass);
            V3 angular_acceleration = T / I;

            R_state deriv;
            deriv.vel = F / s.mass;

            deriv.pos = s.vel;
            
            deriv.ori = s.ori.deriv(s.ang_vel.x, s.ang_vel.y, s.ang_vel.z);
            deriv.ang_vel = angular_acceleration;
            deriv.mass = params.dm;
            deriv.car_mass = params.dm;
            
            return deriv;
        }
        
    };

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


    int main() {
        R_params params;

        R_state s;
        s.pos = V3(0.0f, 0.0f, 0.0f);
        s.vel = V3(0.0f, 0.0f, 0.0f);
        s.ori = quaternion(1.0f, 0.0f, 0.0f, 0.0f);
        s.ang_vel = V3(0, 0, 0);
        s.mass = params.total_mass();
        s.car_mass = params.m_car;
        R_dynamics dynamics(params);
        
        float dt = 0.1f;
        float t = 0.0f;
        //float t_max = 100.0f; 

        for (int i = 0; i < 10000; i++)
        {
            if (i % 1 == 0) {
                float mag_v = s.vel.mag();
                printf("mass : %.2f --- x : %.2f --- y : %.2f --- z : %.2f --- pitch : %.2f ---roll %.2f ---yaw %.2f ---v %.2f ---t : %.2f\n",
                    s.mass, s.pos.x, s.pos.y, s.pos.z, s.pitch / M_PI * 180, s.yaw / M_PI * 180, s.roll / M_PI * 180, mag_v, t);
            }
            
            if (s.mass <= 30.0f && params.Thrust > 0.0f) {
                params.Thrust = 0.0f;
                params.dm = 0.0f;
            }
            

            RK4::rk4(s, dynamics, dt);
            t += dt;
            if (s.pos.y < 0.0f) {
                break;
            }
        }
        return 0;
    };


    //this code contains several mistakes, with bad application of physics principles for inertia and torques and drag , so we need to rework the parametres class and the dynamics class, 
    //this explains why angles do not update correctly and why the rocket does not behave as expected, but it is a good starting point to understand how to structure the code and how to implement the
    //RK4 method for numerical integration of the equations of motion.




