#include "RK4_integrator.hpp"
#include "physics.hpp"
#include <stdio.h>
#include <iostream>
int main() {
    float dt = 0 ;
    float step = 0.01f;              // Time step for the simulation
    float altitude_limit = 0.0f;            // Store the maximum altitude reached
    V3 rotation_angles(0, 0, 0);    // Initialize rotation angles to zero
    rocket_dynamics s;
    s.position = V3(0, 0, 0);
    s.orientation = quaternion(1, 0, 0 , 0); // Example initial orientation
    s.velocity = V3(0, 0, 0);
    s.thrust = V3(0, 0, 1000.0f);   // Example thrust value
    s.mass_fuel = 19.0f;            // Initial fuel mass
    s.angular_velocity = V3(0, 0, 0);
    s.dm = 0.46f;                   // Example mass flow rate

          

    while (s.position.z >= 0.0f && !state) { // Run until the rocket hits the ground or an invalid state is detected
        rotation_angles = s.orientation.euler(); // Update Euler angles for output
        s.drag = s.drag_force(s.velocity, s.position);
        if (s.mass_fuel <= 0.0f) {
            s.thrust = V3(0, 0, 0); // No thrust if fuel is depleted
            s.dm = 0.0f; // Stop fuel consumption
        }
        if (s.position.z > altitude_limit)  {
            altitude_limit = s.position.z; // Update maximum altitude reached
        }
    rk4(s, step);
    std::cout << "step : "<< dt<< "s || tPosition: (" << s.position.x << ", " << s.position.y << ", " << s.position.z << ") \n";
    std::cout << "Velocity: (" << s.velocity.x << ", " << s.velocity.y << ", " << s.velocity.z << ") \n";
    std::cout << "Orientation (euler): (" << rotation_angles.x << ", " << rotation_angles.y << ", " << rotation_angles.z << ") \n";
    std::cout << "thrust: (" << s.thrust.x << ", " << s.thrust.y << ", " << s.thrust.z << " || drag : " << s.drag.x << "," << s.drag.y << "," <<s.drag.z<< " \n";
    std::cout <<"Fuel mass:" << s.mass_fuel<<"kg\n";
    std::cout << "--------------------------------------------------------------" << std::endl;
    dt += step; // Increment time step for the next iteration
}

printf("MAXIMUM REACHED : %.2f m\n", altitude_limit);
    return 0;

}