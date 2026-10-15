#include "quaternion.hpp"
#include <iostream> 
void main_quaternion() {
    // Example usage of the quaternion class
    quaternion q(1, 0.9 , 0.5 , 0.2); // Example quaternion
    q.normalize(); // Normalize the quaternion to ensure it's a valid rotation

    M3 rotation_matrix = q.rotation(); // Get the rotation matrix from the quaternion
    M3 rotation_matrix_transpose = rotation_matrix.transpose(); // Get the transpose of the rotation matrix
    M3 result = rotation_matrix * rotation_matrix_transpose; // Multiply the rotation matrix by
    // Print the rotation matrix
    printf("Identity matrix: %f %f %f\n", result.m[0][0], result.m[1][1], result.m[2][2]);
    
}
