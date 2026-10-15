#pragma once
#include "quaternion.hpp"
#include "rocket_parameters.hpp"
#include "helper_functions.hpp"

class rocket_dynamics {
    public:
    rocket_parameters p;
    quaternion orientation;
    V3 position;
    V3 velocity;
    V3 angular_velocity;
    V3 thrust; // in the body frame
    V3 drag; // in the inertial frame
    float mass_fuel;
    float dm; // fuel mass flow rate
    float get_mass_total(float fuel_mass) const;
    V3 get_CoG(float fuel_mass) const;

    rocket_dynamics();

    M3 get_inertia_tensor(float fuel_mass) const;
    V3 get_d_drag(float fuel_mass) const;
    V3 drag_force(const V3& velocity, const V3& position) const;
    V3 accel(
                const quaternion& orientation,
                const V3& drag_force,
                float fuel_mass) const;
    V3 angular_accel(
                const V3& angular_vel,
                const quaternion& orientation,
                const V3& drag_force,
                float fuel_mass) const;

}; 