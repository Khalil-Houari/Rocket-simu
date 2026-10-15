#pragma once 
#define _USE_MATH_DEFINES
#include "math_tools.hpp"

M3 inertia_cylinder_emp(
    float mass,
    float ext_radius,
    float int_radius,
    float height
);

M3 inertia_cylinder(
    float mass,
    float radius,
    float height
);

M3 inertia_parallelepiped(
    float mass,
    float length,
    float width,
    float height
);

M3 inertia_rhead(
    float density,
    float radius,
    float height
);

M3 steiner(float mass, const V3& d);

float ISA(float z);