#include "cone.h"

Cone::Cone(double r, const Vector & b, const Vector & t)
{
    lower = b;
    upper = t;
    radius = r;
}

double Cone::Radius() const {
    return radius;
}

Vector Cone::Vertex(int n) const {
    return Vector((n & 1) ? upper[0] : lower[0], (n & 2) ? upper[1] : lower[1], (n & 4) ? upper[2] : lower[2]);
}
