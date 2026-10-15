#include "rocket_parameters.hpp"

    // radius_int is a derived parameter, calculated as the external radius minus the thickness of the shell
    float rocket_parameters::radius_int() const { return radius_ext - thickness; }

    // assuming the shape is a perfect ogive , a thing to change;
    V3 rocket_parameters::pos_head() const { return pos_shell + V3(0,0,height); }
    V3 rocket_parameters::CoG_head() const { return V3(0,0,head_height / 3.0f ) + pos_head(); }  

    // assuming the shape is a perfect cylindre , a thing to change
    V3 rocket_parameters::Cog_shell() const { return pos_shell + V3(0,0,height) / 2.0f; } 

    /* the cog is in the middle since the fuel + motor is assumed to be homogeneous,
    but later on when we will rework the motor parameters we will change it*/
    V3 rocket_parameters::CoG_motor() const {return pos_motor + V3(0,0,motor_length) / 2.0f;}
    
    /* for now we will assume that the FC has the same radius as the inner radius of
    the shell, but we can change it later*/
    float rocket_parameters::radius_FC() const { return radius_int(); }  
    
    /* we will assume that the FC is homogeneous*/
    V3 rocket_parameters::CoG_FC() const { return pos_FC + V3(0,0,length_FC) / 2.0f; }
    
    // we will assume that the nanosat is homogeneous and that it is placed at the end of the rocket
    V3 rocket_parameters::CoG_nanosat() const { return pos_nanosat + V3(0,0,height_nanosat) / 2.0f; } 
