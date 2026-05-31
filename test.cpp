#include <iostream>
#include <string>
#include <cmath>
#define M_PI 3.14159265358979323846



struct Point {
float L1 = 2.0f;
float L2 = 1.0f;
float x2,y2;
};

float distance(Point &p, float teta1, float teta2) {
    teta1 = teta1 * M_PI / 180.0f; // Convert to radians
    teta2 = teta2 * M_PI / 180.0f; // Convert to radians
    p.x2 = p.L1 * cos(teta1) + p.L2 * cos(teta2 + teta1);
    p.y2 = p.L1 * sin(teta1) + p.L2 * sin(teta2 + teta1);
    return sqrt(p.x2 * p.x2 + p.y2 * p.y2);
}

int main() {
    Point p;
    float teta1, teta2;
    std::cin >> teta1 >> teta2;
    std::cout << "Distance: " << distance(p, teta1, teta2) << std::endl;
    return 0;
}