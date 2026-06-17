#include <cstdio>
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <algorithm>

float dt; // will be defined , it is not the time step for RK4 
bool state = false; // will be used to break the simulation incase of a bug detection or absurd value 

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

    M3 rotation() const { // rotation matrix
        M3 R;  
        R.m[0][0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
        R.m[0][1] = 2 * (q1 * q2 - q0 * q3);
        R.m[0][2] = 2 * (q1 * q3 + q0 * q2);
        //j
        R.m[1][0] = 2 * (- q1 * q2 + q0 * q3);
        R.m[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
        R.m[1][2] = 2 * (q2 * q3 - q0 * q1);
        //k
        R.m[2][0] = 2 * (q1 * q3 - q0 * q2);
        R.m[2][1] = 2 * (q0 * q1 + q2 * q3);
        R.m[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
        return R;
    }

    V3 euler(V3 angle) { // teta, phi and gamma the rotation angles
        angle.x = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2));
        
        float sin_beta = 2 * (q0 * q2 - q3 * q1);
        angle.y = asin(std::clamp(sin_beta, -0.999999f, 0.999999f)); 

        angle.z = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3));

        return angle;
    }

    quaternion deriv(V3 angular_vel) const {
        float q0Dot = 0.5f * (-q1 * angular_vel.x - q2 * angular_vel.y - q3 * angular_vel.z);
        float q1Dot = 0.5f * (q0 * angular_vel.x + q2 * angular_vel.z - q3 * angular_vel.y);
        float q2Dot = 0.5f * (q0 * angular_vel.y - q1 * angular_vel.z + q3 * angular_vel.x);
        float q3Dot = 0.5f * (q0 * angular_vel.z + q1 * angular_vel.y - q2 * angular_vel.x);
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


// helper functions

M3 inertia_cyclindre_emp(float mass, float ext_radius, float int_radius, float height) {
    M3 inertia;
    inertia.m[0][0] = mass * (3.0f * (ext_radius * ext_radius + int_radius * int_radius) + height * height) / 12.0f; 
    inertia.m[1][1] = inertia.m[0][0]; 
    inertia.m[2][2] =  0.5f * mass * (int_radius * int_radius + ext_radius * ext_radius); // 
    return inertia;
}

M3 inertia_cyclindre(float mass, float radius, float height) {
    M3 inertia;
    inertia.m[0][0] = mass * (3.0f * (radius * radius ) + height * height) / 12.0f; // Ixx = Iyy = (1/12) * m * (3r^2 + h^2)
    inertia.m[1][1] = inertia.m[0][0]; 
    inertia.m[2][2] =  0.5f * mass * (radius * radius ); // Izz = (1/2) * m * r^2
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
    inertia.m[0][0] = M_PI*density*(32.0f * radius * radius * radius * radius * height/315.0f + 29*radius * radius * height * height * height/105.0f); 
    inertia.m[1][1] = inertia.m[0][0];
    inertia.m[2][2] = 0.5f * M_PI * 128.0f * density * radius * radius * radius * radius * height/315.0f; // 
    return inertia;
}
    
M3 steiner(float mass, const V3& d) {
    M3 result;
    
    result.m[0][0] = mass * (d.y * d.y + d.z * d.z); 
    result.m[0][1] = -mass * d.x * d.y;
    result.m[0][2] = -mass * d.x * d.z;

    result.m[1][0] = result.m[0][1]; // Symétrie
    result.m[1][1] = mass * (d.x * d.x + d.z * d.z); 
    result.m[1][2] = -mass * d.y * d.z;

    result.m[2][0] = result.m[0][2]; // Symétrie
    result.m[2][1] = result.m[1][2]; // Symétrie
    result.m[2][2] = mass * (d.x * d.x + d.y * d.y); 

    return result;
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
    V3 pos_shell;
    V3 pos_body;
    V3 pos_head;
    V3 Cog_shell = pos_shell + V3(0,0,height) / 2.0f; // assuming the shape is a perfect cylindre , a thing to change 
    V3 CoG_head = V3(0,0,head_height / 4.0f);
    // motor parameters
    float radius_motor;
    float motor_length;
    float mass_emp_motor;
    float mass_fuel;
    float mass_motor() const {return mass_emp_motor + mass_fuel;}
    float dm; // mass flow rate of the fuel, this one too
    V3 pos_motor;
    V3 CoG_motor() const {return pos_motor + V3(0,0,motor_length) / 2.0f;} // the cog is in the middle since the fuel + motor is assumed to be homogeneous, but later on when we will rework the motor parameters we will change it
    // FC parameters
    float radius_FC = radius_int; // for now we will assume that the FC has the same radius as the inner radius of the shell, but we can change it later 
    float mass_FC;
    float length_FC;
    V3 pos_FC;
    V3 CoG_FC = pos_FC + V3(0,0,length_FC) / 2.0f; // we will assume that the FC is homogeneous
    // nanosat parameters
    float length_nanosat;
    float width_nanosat;
    float height_nanosat;
    float mass_nanosat;
    V3 pos_nanosat;
    V3 CoG_nanosat = pos_nanosat + V3(0,0,length_nanosat) / 2.0f; // we will assume that the nanosat is homogeneous and that it is placed at the end of the rocket
   
    float Cd; // drag coeff
    float Rho; // air density
    float g; // gravity
    V3 CoP;
};

class rocket_state { // still sketching
    public:
    rocket_parameters p;
    quaternion orientation;
    V3 position;
    V3 velocity;

    V3 angular_velocity; 

    V3 thrust; // in the body frame

    V3 drag = V3(0,0,0);   // in the inertial frame
    

    float mass_fuel; // we will update this one during the simulation
    float dm; // mass flow rate of the fuel

    float get_mass_total() const {
        return p.mass_shell + p.mass_motor() + p.mass_FC + p.mass_nanosat;
    }

    V3 get_CoG() const { // this is a function which calculates the CoG with the Huygens Steiner theorem
        return V3( p.Cog_shell * p.mass_shell  + p.CoG_motor() *  p.mass_motor()
               + p.CoG_FC * p.mass_FC  + p.CoG_nanosat * p.mass_nanosat ) / get_mass_total();
    }
 // we will assume that the fuel is distributed uniformly along the motor length
    

    rocket_state() : position(), velocity(), orientation(), angular_velocity(), mass_fuel(0.0f), dm(0.0f) {}

    M3 get_inertia_tensor() const { 
        V3 CoG = get_CoG();
        
        return inertia_cyclindre_emp(p.mass_shell - p.mass_head, p.radius_ext, p.radius_int, p.height-p.head_height)  +steiner(p.mass_shell-p.mass_head,(p.Cog_shell- CoG)) // body shell 

        + inertia_rhead(p.density_shell, p.radius_ext, p.head_height) +steiner(p.mass_head, p.CoG_head- CoG) // head shell

        + inertia_cyclindre(p.mass_FC, p.radius_FC*2, p.length_FC) +steiner(p.mass_FC,(p.CoG_FC-CoG))  // FC

        + inertia_parallelepiped(p.mass_nanosat, p.length_nanosat, p.width_nanosat, p.height_nanosat) +steiner(p.mass_nanosat,(p.CoG_nanosat - CoG)) // nanosat
        
        + inertia_cyclindre(p.mass_motor() ,p.radius_motor,p.motor_length) +steiner(p.mass_motor(), p.CoG_motor() - CoG); // motor
    }

    V3 get_d_drag() const { return p.CoP -  get_CoG();}

    V3 angular_accel(V3 angular_vel) const {
        V3 drag_body =orientation.rotation().transpose() * drag;  // in the body frame      
        V3 d_thr =p.pos_motor  - get_CoG() ;
        V3 d_drag = get_d_drag();
        M3 I = get_inertia_tensor();
        M3 I_inv = I.inverse();
        V3 torque = d_thr.crossV_V(thrust) + d_drag.crossV_V(drag_body);
        return I_inv * (torque - angular_vel.crossV_V(I * angular_vel));
}   
    V3 accel(V3 velocity) const {
        V3 P = V3(0,0,-get_mass_total() * p.g);
        float A = p.radius_ext * p.radius_ext * M_PI;
        V3 drag =  velocity * velocity.magnitude() * p.Rho * p.Cd * A *-0.5f ;
        V3 thrust_inertial = orientation.rotation() * thrust; // in the inertial frame
        V3 accel = ( drag + thrust_inertial + P) / get_mass_total();
        return accel ; 
    }
};

class Rk4 {
    public:
    void rk4_rotational(rocket_state& s, float dt){
        V3 k1 = s.angular_accel(s.angular_velocity);
        V3 k2 = s.angular_accel(s.angular_velocity + k1*dt / 2.0f);
        V3 k3 = s.angular_accel(s.angular_velocity + k2 * dt / 2.0f);
        V3 k4 = s.angular_accel(s.angular_velocity + k3 * dt);
        s.angular_velocity = s.angular_velocity + (k1 + k2 * 2.0f + k3 * 2.0f + k4) * (dt / 6.0f); 
    }

    static void rk4(quaternion& q,rocket_state& s, float dt) {
        quaternion k1 = q.deriv(s.angular_velocity);
        quaternion k2 = (q + k1 * dt / 2.0f).deriv(s.angular_velocity);
        quaternion k3 = (q + k2 * dt / 2.0f).deriv(s.angular_velocity);
        quaternion k4 = (q + k3 * dt).deriv(s.angular_velocity);
        q = q + (k1 + k2 * 2.0f + k3 * 2.0f + k4) * (dt / 6.0f);
        q.normalize();
    }

    void rk4_translational(rocket_state& s, float dt) {
        V3 kv1 = s.accel(s.velocity);
        V3 kv2 = s.accel(s.velocity + kv1 * dt / 2.0f);
        V3 kv3 = s.accel(s.velocity + kv2 * dt / 2.0f);
        V3 kv4 = s.accel(s.velocity + kv3 * dt);
        
        V3 kp1 = s.velocity;
        V3 kp2 = s.velocity + kp1 * dt/2.0f;
        V3 kp3 = s.velocity + kp2 * dt/2.0f;
        V3 kp4 = s.velocity + kp3 * dt;

        s.velocity = s.velocity + (kv1 + kv2 * 2.0f + kv3 * 2.0f + kv4) * (dt / 6.0f);
        s.position = s.position + (kp1 + kp2*2 + kp3*2 + kp4) * (dt/6.0f);
    }
};

int main() {
    for(;;) {

    }
}

    /* 
    *Done :
    -now we have RK4 integrator for the rotation and translational
    -reworked all the algebra tools and some mathematical tools 
    -reworked the quaternion class to be more efficient and more accurate
    -finally finished with inertia tensors and parameters structure even tho we're still missing some parameters 
    -moddeled the the physics
    
    -Priorities :
    *model the atmosphere following the ISA 
    *model the thrust engine
    *execute the loop in the file 


    -notes :
    *we're still treating some variating objects as constants such as the center of pressure 
    *we don't have the drag coefficient graph yet so we assume it to be a constant
    *the model got more accurate since we started not assuming everything to be symmetrical
    *the file will be seperated into many files to ease the maintenance of the code/program
    
    */
   
    /* Last updated : 2026-06-18 */
    



