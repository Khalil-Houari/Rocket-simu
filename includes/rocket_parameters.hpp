#pragma once 
#include "math_tools.hpp"

struct rocket_parameters { // still sketching 
    //shell parameters
    float mass_shell;   // cylindrical shell mass
    float mass_head;    // ogive head mass
    float density_shell;
    float radius_ext;
    float thickness;
    float height;
    float head_height;
    float radius_int() const;
    V3 pos_shell;
    V3 pos_body;
    V3 pos_head() const; // assuming the shape is a perfect ogive , a thing to change;
    V3 Cog_shell() const; // assuming the shape is a perfect cylindre , a thing to change
    V3 CoG_head() const; // assuming the shape is a perfect ogive , a thing to change
    // motor parameters
    float radius_motor;
    float motor_length;
    float mass_emp_motor;
    float mass_fuel;
    float dm = 0.0f; // positive fuel mass flow rate
    V3 pos_motor;
    V3 CoG_motor() const; // the cog is in the middle since the fuel + motor is assumed to be homogeneous, but later on when we will rework the motor parameters we will change it
    // FC parameters
    float radius_FC() const;// for now we will assume that the FC has the same radius as the inner radius of the shell, but we can change it later
    float mass_FC;
    float length_FC;
    V3 pos_FC;
    V3 CoG_FC() const; // we will assume that the FC is homogeneous
    // nanosat parameters
    float length_nanosat;
    float width_nanosat;
    float height_nanosat;
    float mass_nanosat;
    V3 pos_nanosat;
    V3 CoG_nanosat() const; // we will assume that the nanosat is homogeneous and that it is placed at the end of the rocket
   
    float Cd; // drag coeff

    float g; // gravity
    V3 CoP;
};