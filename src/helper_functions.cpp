#include "helper_functions.hpp"
bool state = false;
M3 inertia_cylinder_emp(float mass, float ext_radius, float int_radius, float height) {
    M3 inertia;
    inertia.m[0][0] = mass * (3.0f * (ext_radius * ext_radius + int_radius * int_radius) + height * height) / 12.0f; 
    inertia.m[1][1] = inertia.m[0][0]; 
    inertia.m[2][2] =  0.5f * mass * (int_radius * int_radius + ext_radius * ext_radius); // 
    return inertia;
}

M3 inertia_cylinder(float mass, float radius, float height) {
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


float ISA(float z) {
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