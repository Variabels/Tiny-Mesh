#include "capsule.h"

Capsule::Capsule(double r, double h, const Vector & c)
{
    lower = c;
    upper = lower + Vector(0.0, 0.0, h);
    radius = r;
}

double Capsule::Radius() const {
    return radius;
}

Vector Capsule::Vertex(int n) const {
    return Vector((n & 1) ? upper[0] : lower[0], (n & 2) ? upper[1] : lower[1], (n & 4) ? upper[2] : lower[2]);
}
