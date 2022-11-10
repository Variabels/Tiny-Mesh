#include "cylinder.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

Cylinder::Cylinder(double r, double h, const Vector & c){
    lower = c;
    upper = lower + Vector(h, 0.0, 0.0);
    radius = r;
}

double Cylinder::Radius() const {
    return radius;
}

Vector Cylinder::Vertex(int n) const {
    return Vector((n & 1) ? upper[0] : lower[0], (n & 2) ? upper[1] : lower[1], (n & 4) ? upper[2] : lower[2]);
}
