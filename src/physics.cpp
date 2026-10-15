#include "physics.hpp"

float rocket_dynamics::get_mass_total(float fuel_mass) const {
    return p.mass_shell + p.mass_head + p.mass_emp_motor + std::max(0.0f, fuel_mass)
         + p.mass_FC + p.mass_nanosat;
}

V3 rocket_dynamics::get_CoG(float fuel_mass) const { // this is a function which calculates the CoG with the Huygens Steiner theorem
    float motor_mass = p.mass_emp_motor + std::max(0.0f, fuel_mass);
     return V3( p.Cog_shell() * p.mass_shell  + p.CoG_motor() * motor_mass
         + p.CoG_FC() * p.mass_FC  + p.CoG_nanosat() * p.mass_nanosat + p.CoG_head()* p.mass_head) / get_mass_total(fuel_mass);
}

rocket_dynamics::rocket_dynamics() : p(), orientation(), position(), velocity(), angular_velocity(), thrust(), drag(), mass_fuel(), dm() {
    p.mass_shell = 5.3f;
    p.mass_head = 0.9f;
    p.density_shell = 1800.0f;
    p.radius_ext = 0.10f;
    p.thickness = 0.003f;
    p.height = 1.6f;
    p.head_height = 0.5f;
    p.pos_shell = V3(0, 0, 0);
    p.pos_body = V3(0, 0, 0);

    p.radius_motor = 0.09f;
    p.motor_length = 0.8f;
    p.mass_emp_motor = 2.5f;
    p.mass_fuel = 4.0f;
    p.dm = 0.46f;
    p.pos_motor = V3(0, 0, 0);

    p.mass_FC = 0.4f;
    p.length_FC = 0.15f;
    p.pos_FC = V3(0, 0, 0.9f);

    p.length_nanosat = 0.10f;
    p.width_nanosat = 0.10f;
    p.height_nanosat = 0.10f;
    p.mass_nanosat = 1.33f;
    p.pos_nanosat = V3(0, 0, 1.6f);

    p.Cd = 0.45f;
    p.g = 9.80665f;
    p.CoP = V3(0, 0, 0.45f);

    mass_fuel = p.mass_fuel;
    dm = p.dm;
}

M3 rocket_dynamics::get_inertia_tensor(float fuel_mass) const { 
    float motor_mass = p.mass_emp_motor + std::max(0.0f, fuel_mass);
    V3 CoG = get_CoG(fuel_mass);
    
    return inertia_cylinder_emp(p.mass_shell, p.radius_ext, p.radius_int(), p.height)  +steiner(p.mass_shell,(p.Cog_shell() - CoG)) // body shell 

    + inertia_rhead(p.density_shell, p.radius_ext, p.head_height) +steiner(p.mass_head, p.CoG_head() - CoG) // head shell

    + inertia_cylinder(p.mass_FC, p.radius_FC(), p.length_FC) +steiner(p.mass_FC,(p.CoG_FC() - CoG))  // FC

    + inertia_parallelepiped(p.mass_nanosat, p.length_nanosat, p.width_nanosat, p.height_nanosat) +steiner(p.mass_nanosat,(p.CoG_nanosat() - CoG)) // nanosat
    
    + inertia_cylinder(motor_mass ,p.radius_motor,p.motor_length) +steiner(motor_mass, p.CoG_motor() - CoG); // motor
}

V3 rocket_dynamics::get_d_drag(float fuel_mass) const { return p.CoP - get_CoG(fuel_mass); } // distance between the center of pressure and the center of gravity

V3 rocket_dynamics::drag_force(
    const V3& velocity,
    const V3& position) const {
    constexpr float pi = 3.14159265358979323846f;
    float A = p.radius_ext * p.radius_ext * pi; // cross sectional area
    return velocity * velocity.magnitude() * ISA(position.z) * p.Cd * A * -0.5f;
}

V3 rocket_dynamics::accel(
    const quaternion& orientation,
    const V3& drag_force,
    float fuel_mass) const {
    V3 P = V3(0,0,-get_mass_total(fuel_mass) * p.g);
    V3 drag =  drag_force; // drag force in the inertial frame
    V3 thrust_inertial = orientation.rotation() * thrust; // in the inertial frame
    V3 accel = ( drag + thrust_inertial + P) / get_mass_total(fuel_mass);
    return accel ; 
}

V3 rocket_dynamics::angular_accel(
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

