#define _USE_MATH_DEFINES
#include <cstdio>
#include <iostream>
#include <cmath>
#include <algorithm>


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
// thu
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
        R.m[1][0] = 2 * (q1 * q2 + q0 * q3);
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
    float mass_shell;   // cylindrical shell mass
    float mass_head;    // ogive head mass
    float density_shell;
    float radius_ext;
    float thickness;
    float height;
    float head_height;
    float radius_int() const { return radius_ext - thickness; }
    V3 pos_shell;
    V3 pos_body;
    V3 pos_head() const { return pos_shell + V3(0,0,height); } // assuming the shape is a perfect ogive , a thing to change;
    V3 Cog_shell() const { return pos_shell + V3(0,0,height) / 2.0f; } // assuming the shape is a perfect cylindre , a thing to change
    V3 CoG_head() const { return V3(0,0,head_height / 3.0f ) + pos_head(); } // assuming the shape is a perfect ogive , a thing to change
    // motor parameters
    float radius_motor;
    float motor_length;
    float mass_emp_motor;
    float mass_fuel;
    float dm = 0.0f; // positive fuel mass flow rate
    V3 pos_motor;
    V3 CoG_motor() const {return pos_motor + V3(0,0,motor_length) / 2.0f;} // the cog is in the middle since the fuel + motor is assumed to be homogeneous, but later on when we will rework the motor parameters we will change it
    // FC parameters
    float radius_FC() const { return radius_int(); } // for now we will assume that the FC has the same radius as the inner radius of the shell, but we can change it later
    float mass_FC;
    float length_FC;
    V3 pos_FC;
    V3 CoG_FC() const { return pos_FC + V3(0,0,length_FC) / 2.0f; } // we will assume that the FC is homogeneous
    // nanosat parameters
    float length_nanosat;
    float width_nanosat;
    float height_nanosat;
    float mass_nanosat;
    V3 pos_nanosat;
    V3 CoG_nanosat() const { return pos_nanosat + V3(0,0,height_nanosat) / 2.0f; } // we will assume that the nanosat is homogeneous and that it is placed at the end of the rocket
   
    float Cd; // drag coeff

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


    float ISA(float z) const {
        float g = 9.80665f; // standard gravity
        float R = 287.05f; // specific gas constant for dry air
        float r = 6356766.0f; // radius of the earth
        float hg = (r*z)/(r+z); // geopotential altitude
        float d, T, P;
        //troposphere
        float hb = 0.0f;
        float Tb = 288.15f;
        float Lb = -0.0065f;
        float Pb = 101325.0f;
        //tropopause
        if (hg >= 11000.0f && hg < 20000.0f) { 
            hb = 11000.0f;
            Tb = 216.65f;
            Lb = 0.000f;
            Pb = 22632.1f;
        }
        //stratosphere (lower)
        else if (hg >= 20000.0f && hg < 32000.0f) { 
            hb = 20000.0f;
            Tb = 216.65f;
            Lb = 0.0010f;
            Pb = 5474.89f;
        }
        //stratosphere (upper)
        else if (hg >= 32000.0f && hg < 47000.0f) { 
            hb = 32000.0f;
            Tb = 228.65f;
            Lb = 0.0028f;
            Pb = 868.02f;
        }
        //stratopause
        else if (hg >= 47000.0f && hg < 51000.0f){ 
            hb = 47000.0f;
            Tb = 270.65f;
            Lb = 0.0f;
            Pb = 110.91f;
        }
        // isothermal layer
        if (Lb == 0.0f) {
            T = Tb;
            P = Pb * exp(-1.0f * g * (hg- hb)/ (Tb * R));
        } else {
        // gradient layer
        T = Tb + Lb * (hg - hb);
        P = Pb * pow(Tb / T, g / (Lb * R));
    }
    d = P / (R * T); // density
    return d; // Note : we can notice that we're only returning the density, 
            //but we can also return the temperature and pressure if needed
    }

    float mass_fuel; // we will update this one during the simulation
    float dm; // mass flow rate of the fuel

    float get_mass_total(float fuel_mass) const {
        return p.mass_shell + p.mass_head + p.mass_emp_motor + std::max(0.0f, fuel_mass)
             + p.mass_FC + p.mass_nanosat;
    }

    V3 get_CoG(float fuel_mass) const { // this is a function which calculates the CoG with the Huygens Steiner theorem
        float motor_mass = p.mass_emp_motor + std::max(0.0f, fuel_mass);
         return V3( p.Cog_shell() * p.mass_shell  + p.CoG_motor() * motor_mass
             + p.CoG_FC() * p.mass_FC  + p.CoG_nanosat() * p.mass_nanosat + p.CoG_head()* p.mass_head) / get_mass_total(fuel_mass);
    }
 // we will assume that the fuel is distributed uniformly along the motor length
    

    rocket_state() : p(), orientation(), position(), velocity(), angular_velocity(), mass_fuel(0.0f), dm(0.0f), drag() {
        // coque (cylindre) et ogive : composite, épaisseur 3 mm
        p.mass_shell = 5.3f;
        p.mass_head = 0.9f;
        p.density_shell = 1800.0f;
        p.radius_ext = 0.10f;
        p.thickness = 0.003f;
        p.height = 1.6f;
        p.head_height = 0.5f;
        p.pos_shell = V3(0, 0, 0);

        // moteur (en bas)
        p.radius_motor = 0.09f;
        p.motor_length = 0.8f;
        p.mass_emp_motor = 2.5f;
        p.mass_fuel = 4.0f;
        p.dm = 0.46f;
        p.pos_motor = V3(0, 0, 0);

        // calculateur de vol (au-dessus du moteur)
        p.mass_FC = 0.4f;
        p.length_FC = 0.15f;
        p.pos_FC = V3(0, 0, 0.9f);

        // nanosat 1U, à la base de l'ogive
        p.length_nanosat = 0.10f;
        p.width_nanosat = 0.10f;
        p.height_nanosat = 0.10f;
        p.mass_nanosat = 1.33f;
        p.pos_nanosat = V3(0, 0, 1.6f);

        p.Cd = 0.45f;
        p.g = 9.80665f;
        p.CoP = V3(0, 0, 0.45f);
    }

    M3 get_inertia_tensor(float fuel_mass) const { 
        float motor_mass = p.mass_emp_motor + std::max(0.0f, fuel_mass);
        V3 CoG = get_CoG(fuel_mass);
        
        return inertia_cyclindre_emp(p.mass_shell, p.radius_ext, p.radius_int(), p.height)  +steiner(p.mass_shell,(p.Cog_shell() - CoG)) // body shell 

        + inertia_rhead(p.density_shell, p.radius_ext, p.head_height) +steiner(p.mass_head, p.CoG_head() - CoG) // head shell

        + inertia_cyclindre(p.mass_FC, p.radius_FC(), p.length_FC) +steiner(p.mass_FC,(p.CoG_FC() - CoG))  // FC

        + inertia_parallelepiped(p.mass_nanosat, p.length_nanosat, p.width_nanosat, p.height_nanosat) +steiner(p.mass_nanosat,(p.CoG_nanosat() - CoG)) // nanosat
        
        + inertia_cyclindre(motor_mass ,p.radius_motor,p.motor_length) +steiner(motor_mass, p.CoG_motor() - CoG); // motor
    }

    // helper functions for drag 
    V3 drag_force(
        const V3& velocity,
        const V3& position) const {
        float A = p.radius_ext * p.radius_ext * M_PI; // cross sectional area
        return velocity * velocity.magnitude() * ISA(position.z) * p.Cd * A * -0.5f;
    }
    V3 get_d_drag(float fuel_mass) const { return p.CoP - get_CoG(fuel_mass); } // distance between the center of pressure and the center of gravity


    V3 angular_accel(
        const V3& angular_vel, 
        const quaternion& orientation, 
        const V3& drag_force,
        float fuel_mass) const {
        V3 drag_body =orientation.rotation().transpose() * drag_force;  // in the body frame      
        V3 d_thr = p.pos_motor - get_CoG(fuel_mass);
        V3 d_drag = get_d_drag(fuel_mass);
        M3 I = get_inertia_tensor(fuel_mass);
        M3 I_inv = I.inverse();
        V3 torque = d_thr.crossV_V(thrust) + d_drag.crossV_V(drag_body);
        return I_inv * (torque - angular_vel.crossV_V(I * angular_vel));
}   
    V3 accel(
        const quaternion& orientation,
        const V3& drag_force,
        float fuel_mass) const {
        V3 P = V3(0,0,-get_mass_total(fuel_mass) * p.g);
        V3 drag =  drag_force; // drag force in the inertial frame
        V3 thrust_inertial = orientation.rotation() * thrust; // in the inertial frame
        V3 accel = ( drag + thrust_inertial + P) / get_mass_total(fuel_mass);
        return accel ; 
    }
};

// Runge-Kutta 4th order integration for the rocket state
    void rk4(rocket_state& s, float dt){
        
        // initialization of the RK4 stages
        float fuel_mass_init = s.mass_fuel;
        const float fuel_dm = (s.dm > 0.0f) ? -s.dm : s.dm; // If dm is stored as a positive burn magnitude, make it negative so fuel decreases.
        V3 w0 = s.angular_velocity;
        quaternion q0 = s.orientation;
        V3 init_stage_p = s.position;
        V3 init_stage_vel = s.velocity;
        
        
        // Stage 1
        float fuel_mass_stage = fuel_mass_init;
        quaternion stage_q = q0;
        V3 stage_w = w0;
        V3 stage = init_stage_p;
        V3 stage_vel = init_stage_vel;
        V3 stage_drag = s.drag_force(stage_vel, stage);

        float km1 = fuel_dm;
        V3 Kv1 = s.accel(stage_q, stage_drag, fuel_mass_stage);
        V3 Kp1 = stage_vel;
        quaternion kq1 = stage_q.deriv(stage_w);
        V3 kw1 = s.angular_accel(stage_w, stage_q, stage_drag, fuel_mass_stage);

        
        //stage 2
        fuel_mass_stage = fuel_mass_init + km1 * dt / 2.0f;
        stage_q = q0 + kq1 * dt / 2.0f;
        stage_q.normalize(); // Normalize the quaternion to prevent drift
        stage_w = w0 + kw1 * dt / 2.0f;
        stage = init_stage_p + Kp1 * dt / 2.0f;
        stage_vel = init_stage_vel + Kv1 * dt / 2.0f;
        stage_drag = s.drag_force(stage_vel, stage);

        float km2 = fuel_dm;
        V3 Kv2 = s.accel(stage_q, stage_drag, fuel_mass_stage);
        V3 Kp2 = stage_vel;
        quaternion kq2 = stage_q.deriv(stage_w);
        V3 kw2 = s.angular_accel(stage_w, stage_q, stage_drag, fuel_mass_stage);
        

        //stage 3
        fuel_mass_stage = fuel_mass_init + 0.5f * km2 * dt;
        stage_q = q0 + kq2 * dt / 2.0f;
        stage_q.normalize(); // Normalize the quaternion to prevent drift
        stage_w = w0 + kw2 * dt / 2.0f;
        stage = init_stage_p + Kp2 * dt / 2.0f;
        stage_vel = init_stage_vel + Kv2 * dt / 2.0f;
        stage_drag = s.drag_force(stage_vel, stage);
        
        float km3 = fuel_dm;
        V3 Kv3 = s.accel(stage_q, stage_drag, fuel_mass_stage);
        V3 Kp3 = stage_vel;
        quaternion kq3 = stage_q.deriv(stage_w);
        V3 kw3 = s.angular_accel(stage_w, stage_q, stage_drag, fuel_mass_stage);
        

        //stage 4
        fuel_mass_stage = fuel_mass_init + km3 * dt;
        stage_q = q0 + kq3 * dt;
        stage_q.normalize(); // Normalize the quaternion to prevent drift
        stage_w = w0 + kw3 * dt;
        stage = init_stage_p + Kp3 * dt;
        stage_vel = init_stage_vel + Kv3 * dt;
        stage_drag = s.drag_force(stage_vel, stage);

        float km4 = fuel_dm;
        V3 Kv4 = s.accel(stage_q, stage_drag, fuel_mass_stage);
        V3 Kp4 = stage_vel;
        quaternion kq4 = stage_q.deriv(stage_w);
        V3 kw4 = s.angular_accel(stage_w, stage_q, stage_drag, fuel_mass_stage);
        

        // Final update
        s.mass_fuel = fuel_mass_init + (km1 + 2.0f * km2 + 2.0f * km3 + km4) * (dt / 6.0f);
        s.orientation = s.orientation+ (kq1 + kq2 * 2.0f + kq3 * 2.0f + kq4) * (dt / 6.0f);
        s.angular_velocity = s.angular_velocity + (kw1 + kw2 * 2.0f + kw3 * 2.0f + kw4) * (dt / 6.0f); 
        s.velocity = s.velocity + (Kv1 + Kv2*2.0 + Kv3*2.0 + Kv4) * (dt / 6.0f);
        s.position = s.position + (Kp1 + Kp2*2 + Kp3*2 + Kp4) * (dt/6.0f);
        s.orientation.normalize();
    }

int main() {
    float dt = 0 ;
    float step = 1.0f; // Time step for the simulation
    V3 rotation_angles(0, 0, 0); // Initialize rotation angles to zero
    rocket_state s;
    s.position = V3(0, 0, 0);
    s.orientation = quaternion(1, 0, 0, 0);
    s.velocity = V3(0, 0, 0);
    s.thrust = V3(0, 0, 1000.0f); // Example thrust value
    s.mass_fuel = 19.0f; // Initial fuel mass
    s.angular_velocity = V3(0, 0, 0);
    s.dm = 0.46f; // Example mass flow rate

    M3 rot = s.orientation.rotation();
    M3 rot_transpose = rot.transpose();
    M3 result = rot * rot_transpose;

    printf(result.m[0][0] == 1.0f && result.m[1][1] == 1.0f && result.m[2][2] == 1.0f ? "Rotation matrix is orthogonal.\n" : "Rotation matrix is NOT orthogonal.\n");
            /* /////////////////////////////

    while (s.position.z >= 0.0f && !state) { // Run until the rocket hits the ground or an invalid state is detected
        s.orientation.euler(rotation_angles); // Update Euler angles for output
        s.drag = s.drag_force(s.velocity, s.position);
        if (s.mass_fuel <= 0.0f) {
            s.thrust = V3(0, 0, 0); // No thrust if fuel is depleted
            s.dm = 0.0f; // Stop fuel consumption
        }
    rk4(s, step);
    std::cout << "step : "<< dt<< "s || tPosition: (" << s.position.x << ", " << s.position.y << ", " << s.position.z << ") \n";
    std::cout << "Velocity: (" << s.velocity.x << ", " << s.velocity.y << ", " << s.velocity.z << ") \n";
    std::cout << "Orientation (euler): (" << rotation_angles.x << ", " << rotation_angles.y << ", " << rotation_angles.z << ") \n";
    std::cout << "thrust: (" << s.thrust.x << ", " << s.thrust.y << ", " << s.thrust.z << " || drag : " << s.drag.x << "," << s.drag.y << "," <<s.drag.z<< " \n";
    std::cout <<"Fuel mass:" << s.mass_fuel<<"kg\n";
    std::cout << "--------------------------------------------------------------" << std::endl;
    dt += 0.1f; // Increment time step for the next iteration
}

            /////////////////////*/ 
}

/*
 * Done:
         - V3, M3 and quaternion implemented.
         - Quaternion rotation matrix, derivative and Euler-angle conversion implemented.
         - Rocket parameters now use const getter methods for derived geometry and CoGs.
         - Ogive shell, internal payload/mousse, flight computer, motor and fuel are
             represented as separate mass contributors in the CoG calculation.
         - ISA density model implemented through the 47-51 km layer using geopotential altitude.
         - Still-air drag law implemented.
         - RK4 couples position, velocity, orientation, angular velocity and fuel mass.
         - Constant fuel mass-flow rate is included in the RK4 stages.
         - Orientation has one source of truth: s.orientation.
         - A terminal simulation loop and basic state output exist.
         - Test R * R^T = I for the rotation matrix .

 * Missing physics:
         - Wind model and relative air velocity.
         - Angle of attack from body-axis and relative-air-velocity directions.
         - Lift, normal force, aerodynamic side force and aerodynamic moments.
         - Mach/Reynolds-dependent aerodynamic coefficients and a Cd lookup/model.
         - Pitch/yaw damping and angle-of-attack-dependent center of pressure.
         - Engine model: thrust curve, mass-flow curve, burn duration, shutdown and thrust direction.
         - Remaining ISA layers above 51 km, if higher-altitude flight is required but not urgent.
         - Consistent ogive-shell inertia and internal payload/mousse inertia validation.
         - dI/dt effects during propellant burn.

 * Required corrections before trusting the simulation:
         - Initialize s.mass_fuel and s.dm from the rocket parameters.
         - Clamp fuel mass and burn rate when the remaining fuel is less than dm * dt.
         - Refresh s.drag after every integration step if it is kept as stored state/output.
         - Make inverse() return safely when the inertia matrix is singular.
         - Guard quaternion normalization against a zero or non-finite norm.
         - Add finite-value, positive-mass and valid-inertia state checks after every step.
         - Add ground-contact handling so position.z cannot remain below zero.

 * Validation:
         - Test rotation matrices with identity and known 90-degree rotations.
         - Test free fall with thrust and drag disabled.
         - Test constant-thrust vertical flight.
         - Test constant-burn fuel depletion against m(t) = m0 - dm * t.
         - Test conservation of energy/angular momentum where applicable.

 * Later:
         - Separate math, rocket physics, integrator and main into different files.
         - Export trajectory data to CSV and plot it before adding a GUI/OpenGL layer.
*/

/* Last updated: 2026-10-14 */
   

    



